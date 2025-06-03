#!/usr/bin/env python3
import subprocess
import os
import json
import logging
import tempfile
from pathlib import Path
from typing import Dict, Any, List

from utils import setup_directories, validate_input, setup_logging


class SlurmPipelineRSV:
    def __init__(self, config_path: str = "config/paths.json"):
        self.config = self.load_config(config_path)
        self.references = self.config["references"]
        self.directories = self.config["data"]["output"]
        self.slurm_params = self.config["slurm"]
        self.modules = self.config["modules"]

        # Création des dossiers de sortie
        setup_directories(self.directories)

        # Configuration du logging
        setup_logging(self.config["logs"])
        logging.info("Initialisation du pipeline RSV avec Slurm")

    def load_config(self, config_path: str) -> Dict[str, Any]:
        """Charge la configuration depuis le fichier JSON"""
        with open(config_path) as f:
            return json.load(f)

    def create_slurm_script(self, commands: List[str], job_name: str, dependencies: List[str] = None) -> str:
        """Crée un script Slurm temporaire"""
        # Générer les lignes pour charger les modules
        module_lines = "\n".join(f"module load {mod}" for mod in self.modules.values())
        
        # Ligne pour les dépendances (si présentes)
        dependency_line = ""
        if dependencies:
            dependency_line = f"#SBATCH --dependency=afterok:{':'.join(dependencies)}"

        # Préparer les commandes avec des sauts de ligne AVANT le f-string
        commands_str = "\n".join(commands)

        # Créer le template final sans erreur de syntaxe
        slurm_template = f"""#!/bin/bash
    #SBATCH --account={self.slurm_params["account"]}
    #SBATCH --job-name={job_name}
    #SBATCH --output={self.config["logs"]}/{job_name}_%j.out
    #SBATCH --error={self.config["logs"]}/{job_name}_%j.err
    #SBATCH --time={self.slurm_params["time"]}
    #SBATCH --mem={self.slurm_params["mem"]}
    #SBATCH --cpus-per-task={self.slurm_params["cpus"]}
    #SBATCH --partition={self.slurm_params["partition"]}
    {dependency_line}
    # Chargement des modules
    module purge
    {module_lines}

    # Commandes
    {commands_str}
    """

        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(slurm_template)
            return f.name

    def submit_slurm_job(self, commands: List[str], step_name: str, dependencies: List[str] = None) -> str:
        """Soumet un job Slurm et retourne l'ID du job"""
        try:
            script_path = self.create_slurm_script(commands, step_name, dependencies)
            cmd = f"sbatch --parsable {script_path}"
            logging.info(f"Soumission Slurm ({step_name}): {cmd}")
            result = subprocess.check_output(cmd, shell=True, text=True).strip()
            os.unlink(script_path)
            return result
        except subprocess.CalledProcessError as e:
            logging.error(f"Échec soumission Slurm ({step_name}): {e}")
            raise

    def concatenate_fastq(self) -> List[str]:
        jobs = []
        input_dir = Path(self.config["data"]["input"])
        output_dir = Path(self.directories["concatenated"])
        output_dir.mkdir(parents=True, exist_ok=True)
        for barcode_dir in input_dir.iterdir():
            if barcode_dir.is_dir():
                output_file = output_dir / f"{barcode_dir.name}.fastq.gz"
                if output_file.exists():
                    logging.info(f"Fichier {output_file} existe déjà - skip")
                    continue
                cmd = f"zcat {barcode_dir}/*.fastq.gz | gzip > {output_file}"
                job_id = self.submit_slurm_job([cmd], f"concat_{barcode_dir.name}")
                jobs.append(job_id)
        return jobs

    def trim_adapters(self, concat_job_ids: List[str]) -> List[str]:
        jobs = []
        input_dir = Path(self.directories["concatenated"])
        output_dir = Path(self.directories["trimmed"])
        output_dir.mkdir(parents=True, exist_ok=True)
        input_parent = Path(self.config["data"]["input"])
        expected_files = [input_dir / f"{d.name}.fastq.gz" for d in input_parent.iterdir() if d.is_dir()]
        for fq in expected_files:
            base = fq.stem
            output_file = output_dir / f"{base}_trimmed.fastq.gz"
            if output_file.exists():
                logging.info(f"Fichier trimmed {output_file} existe déjà - skip")
                continue
            cmd = f"cutadapt -u 20 -u -20 -o {output_file} {fq}"
            job_id = self.submit_slurm_job([cmd], f"cutadapt_{base}", dependencies=concat_job_ids)
            jobs.append(job_id)
        return jobs

    def run_nanoplot(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_dir = Path(self.directories["trimmed"])
        output_dir = Path(self.directories["nanoplot"])
        output_dir.mkdir(parents=True, exist_ok=True)
        for fq in input_dir.glob('*_trimmed.fastq.gz'):
            base = fq.name.replace('_trimmed.fastq.gz', '')
            outdir = output_dir / f"{base}_quality"
            if (outdir / "NanoPlot-report.html").exists():
                logging.info(f"NanoPlot pour {base} déjà fait - skip")
                continue
            cmd = f"NanoPlot --fastq {fq} --outdir {outdir}"
            job_id = self.submit_slurm_job([cmd], f"nanoplot_{base}", dependencies=dependencies)
            jobs.append(job_id)
        return jobs

    def run_minimap2(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_dir = Path(self.directories["trimmed"])
        output_dir = Path(self.directories["aligned"])
        output_dir.mkdir(parents=True, exist_ok=True)

        for fq in input_dir.glob('*_trimmed.fastq.gz'):
            base = fq.name.replace('_trimmed.fastq.gz', '')
            sam_a = output_dir / f"{base}_refA.sam"
            sam_b = output_dir / f"{base}_refB.sam"

            # Vérifier si les fichiers SAM existent déjà
            if sam_a.exists() and sam_b.exists():
                logging.info(f"Alignements pour {base} existent déjà → skip")
                continue

            # Lancer seulement les alignements manquants
            cmd_a = f"minimap2 -ax map-ont {self.references['rsv_a']} {fq} > {sam_a}"
            cmd_b = f"minimap2 -ax map-ont {self.references['rsv_b']} {fq} > {sam_b}"

            job_a = self.submit_slurm_job([cmd_a], f"minimapA_{base}", dependencies=dependencies)
            job_b = self.submit_slurm_job([cmd_b], f"minimapB_{base}", dependencies=dependencies)
            jobs.extend([job_a, job_b])

        return jobs

    def process_sam_files(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_dir = Path(self.directories["aligned"])
        output_dir = Path(self.directories["sorted_bam"])
        output_dir.mkdir(parents=True, exist_ok=True)
        for sam in input_dir.glob('*.sam'):
            base = sam.stem
            bam = output_dir / f"{base}_sorted.bam"
            if bam.exists():
                logging.info(f"BAM {bam} existe déjà - skip")
                continue
            cmds = [
                f"samtools view -b {sam} | samtools sort -o {bam}",
                f"samtools index {bam}"
            ]
            job_id = self.submit_slurm_job(cmds, f"samtools_{base}", dependencies=dependencies)
            jobs.append(job_id)
        return jobs

    def run_genome_coverage(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_dir = Path(self.directories["sorted_bam"])
        output_dir = Path(self.directories["genomecov"])
        output_dir.mkdir(parents=True, exist_ok=True)
        for bam in input_dir.glob('*_sorted.bam'):
            bed = output_dir / f"{bam.stem}_genomecov.bed"
            if bed.exists():
                logging.info(f"Couverture {bed} existe déjà - skip")
                continue
            cmd = f"bedtools genomecov -ibam {bam} -bga > {bed}"
            job_id = self.submit_slurm_job([cmd], f"gencov_{bam.stem}", dependencies=dependencies)
            jobs.append(job_id)
        return jobs

    def extract_low_coverage(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_dir = Path(self.directories["genomecov"])
        output_dir = Path(self.directories["extractline"])
        output_dir.mkdir(parents=True, exist_ok=True)
        for bed in input_dir.glob('*_genomecov.bed'):
            base = bed.stem.replace('_genomecov', '')
            out_bed = output_dir / f"{base}_low.bed"
            if out_bed.exists():
                logging.info(f"Fichier low coverage {out_bed} existe déjà - skip")
                continue
            cmd = f"awk '$4 <= 30 {{print}}' {bed} > {out_bed}"
            job_id = self.submit_slurm_job([cmd], f"extract_{base}", dependencies=dependencies)
            jobs.append(job_id)
        return jobs

    def call_variants(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_dir = Path(self.directories["sorted_bam"])
        output_dir = Path(self.directories["medaka_variants"])
        output_dir.mkdir(parents=True, exist_ok=True)
        for bam in input_dir.glob('*_sorted.bam'):
            base = bam.stem.replace('_sorted', '')
            ref_key = 'rsv_a' if 'refA' in base else 'rsv_b'
            out_dir = output_dir / f"{base}_medaka_variant"
            if (out_dir / "medaka.sorted.vcf").exists():
                logging.info(f"Variants {out_dir} existants - skip")
                continue
            cmd = f"medaka_variant -i {bam} -r {self.references[ref_key]} -o {out_dir}"
            job_id = self.submit_slurm_job([cmd], f"medaka_{base}", dependencies=dependencies)
            jobs.append(job_id)
        return jobs

    def compress_and_index_vcf(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_dir = Path(self.directories["medaka_variants"])
        for vcf in input_dir.rglob('medaka.sorted.vcf'):
            vcf_gz = vcf.with_suffix('.vcf.gz')
            tbi_file = vcf.with_suffix('.vcf.gz.tbi')

            # Vérifie si les fichiers existent déjà
            if vcf_gz.exists() and tbi_file.exists():
                logging.info(f"Fichiers {vcf_gz} et {tbi_file} existent déjà → skip")
                continue

            # Commandes à exécuter
            cmds = [
                f"bgzip -f {vcf}",           # -f force overwrite si fichier existe
                f"bcftools index {vcf}.gz"   # Génère l'index .tbi
            ]

            # Soumission du job Slurm
            job_id = self.submit_slurm_job(cmds, f"compress_{vcf.parent.name}", dependencies=dependencies)
            jobs.append(job_id)

        return jobs

    def generate_consensus(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_dir = Path(self.directories["medaka_variants"])
        extract_dir = Path(self.directories["extractline"])
        consensus_dir = Path(self.directories["consensus"])
        consensus_dir.mkdir(parents=True, exist_ok=True)

        for vcf in input_dir.rglob('medaka.sorted.vcf.gz'):
            # Extraire proprement le préfixe
            folder_name = vcf.parent.name  # ex: barcode14.fastq_refA_medaka_variant
            base_part = folder_name.split('_medaka_variant')[0]  # ex: barcode14.fastq_refA
            ref_type = 'refA' if '_refA' in base_part else 'refB'
            base = base_part.split('.')[0]  # ex: barcode14

            # Construire le chemin vers le fichier BED attendu
            bed_file = extract_dir / f"{base_part}_sorted_low.bed"  # ex: barcode14.fastq_refA_sorted_low.bed

            # Fichier de sortie
            ref_key = 'rsv_a' if ref_type == 'refA' else 'rsv_b'
            consensus_file = consensus_dir / f"{base}_consensus_{ref_type}.fasta"

            if consensus_file.exists():
                logging.info(f"Consensus {consensus_file} existe déjà - skip")
                continue

            if not bed_file.exists():
                logging.warning(f"Fichier BED manquant pour {base} → consensus non lancé")
                continue

            # Commande bcftools
            cmd = (
                f"bcftools consensus -m {bed_file} "
                f"-f {self.references[ref_key]} -o {consensus_file} {vcf}"
            )

            job_id = self.submit_slurm_job([cmd], f"consensus_{base}_{ref_type}", dependencies=dependencies)
            jobs.append(job_id)

        return jobs

    def execute_pipeline(self):
        """Exécute toutes les étapes du pipeline en respectant les dépendances"""
        try:
            logging.info("Début du pipeline")
            concat_jobs = self.concatenate_fastq()
            trim_jobs = self.trim_adapters(concat_jobs)
            nanoplot_jobs = self.run_nanoplot(trim_jobs)
            minimap_jobs = self.run_minimap2(trim_jobs)
            sam_jobs = self.process_sam_files(minimap_jobs)
            cov_jobs = self.run_genome_coverage(sam_jobs)
            extract_jobs = self.extract_low_coverage(cov_jobs)
            medaka_jobs = self.call_variants(extract_jobs)
            compress_jobs = self.compress_and_index_vcf(medaka_jobs)
            consensus_jobs = self.generate_consensus(compress_jobs)
            logging.info(f"Pipeline terminé. Consensus généré dans {self.directories['consensus']}")
        except Exception as e:
            logging.critical(f"Échec du pipeline: {str(e)}")
            raise


if __name__ == "__main__":
    pipeline = SlurmPipelineRSV()
    pipeline.execute_pipeline()