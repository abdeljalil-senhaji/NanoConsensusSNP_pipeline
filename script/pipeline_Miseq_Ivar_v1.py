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
    def __init__(self, config_path: str = "config/paths2.json"):
        self.config = self.load_config(config_path)
        self.references = self.config["references"]
        self.directories = self.config["data"]["output"]
        self.slurm_params = self.config["slurm"]
        self.modules = self.config["modules"]

        # Création des dossiers de sortie
        setup_directories(self.directories)

        # Configuration du logging
        setup_logging(self.config["logs"])

        logging.info("Initialisation du pipeline RSV Illumina avec Slurm")

    def load_config(self, config_path: str) -> Dict[str, Any]:
        """Charge la configuration depuis le fichier JSON"""
        with open(config_path) as f:
            return json.load(f)

    def create_slurm_script(self, commands: List[str], job_name: str, dependencies: List[str] = None) -> str:
        """Crée un script Slurm temporaire"""

        module_lines = "\n".join(f"module load {mod}" for mod in self.modules.values())
        dependency_line = ""
        if dependencies:
            dependency_line = f"#SBATCH --dependency=afterok:{':'.join(dependencies)}"

        # Préparer les commandes avec des sauts de ligne AVANT le f-string
        commands_str = "\n".join(commands)

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

    def trim_adapters(self) -> List[str]:
        """Trim les adaptateurs Illumina sur R1 et R2"""
        jobs = []
        input_dir = Path(self.config["data"]["input"])
        output_dir = Path(self.directories["trimmed"])
        output_dir.mkdir(parents=True, exist_ok=True)
    
        # Recherche tous les fichiers R1 et R2
        r1_files = {f.name.split('_R1')[0]: f for f in input_dir.glob('*_R1*.fastq.gz')}
        r2_files = {f.name.split('_R2')[0]: f for f in input_dir.glob('*_R2*.fastq.gz')}
    
        for base, r1 in r1_files.items():
            r2 = r2_files.get(base)
            if not r2:
                logging.warning(f"Pas de fichier R2 trouvé pour {base} → skip")
                continue
    
            out_r1 = output_dir / f"{base}_R1_trimmed.fastq.gz"
            out_r2 = output_dir / f"{base}_R2_trimmed.fastq.gz"
    
            if out_r1.exists() and out_r2.exists():
                logging.info(f"Fichiers trimmed pour {base} déjà présents → skip")
                continue
    
            cmd = f"cutadapt -u 20 -u -20 -o {out_r1} -p {out_r2} {r1} {r2}"
            job_id = self.submit_slurm_job([cmd], f"cutadapt_{base}")
            jobs.append(job_id)
    
        return jobs

    def run_bowtie2(self, dependencies: List[str]) -> List[str]:
        """Alignement avec bowtie2 sur RSV A et RSV B"""
        jobs = []
        input_dir = Path(self.directories["trimmed"])
        output_dir = Path(self.directories["aligned"])
        output_dir.mkdir(parents=True, exist_ok=True)

        for r1 in input_dir.glob('*_R1_trimmed.fastq.gz'):
            base = r1.name.replace('_R1_trimmed.fastq.gz', '')
            r2 = input_dir / f"{base}_R2_trimmed.fastq.gz"
            if not r2.exists():
                logging.warning(f"Fichier R2 manquant pour {base} → skip")
                continue

            sam_a = output_dir / f"{base}_refA.sam"
            sam_b = output_dir / f"{base}_refB.sam"

            if sam_a.exists() and sam_b.exists():
                logging.info(f"Alignements pour {base} existent déjà → skip")
                continue

            cmd_a = f"bowtie2 -x {self.references['rsv_a_index']} -1 {r1} -2 {r2} -S {sam_a}"
            cmd_b = f"bowtie2 -x {self.references['rsv_b_index']} -1 {r1} -2 {r2} -S {sam_b}"

            job_a = self.submit_slurm_job([cmd_a], f"bt2_A_{base}", dependencies=dependencies)
            job_b = self.submit_slurm_job([cmd_b], f"bt2_B_{base}", dependencies=dependencies)

            jobs.extend([job_a, job_b])

        return jobs

    def process_sam_files(self, dependencies: List[str]) -> List[str]:
        """Conversion SAM → BAM trié + indexation"""
        jobs = []
        input_dir = Path(self.directories["aligned"])
        output_dir = Path(self.directories["sorted_bam"])
        output_dir.mkdir(parents=True, exist_ok=True)

        for sam in input_dir.glob('*.sam'):
            base = sam.stem
            bam = output_dir / f"{base}.bam"
            sorted_bam = output_dir / f"{base}_sorted.bam"

            if sorted_bam.exists():
                logging.info(f"BAM trié {sorted_bam} existe déjà → skip")
                continue

            cmds = [
                f"samtools view -b {sam} > {bam}",
                f"samtools sort {bam} -o {sorted_bam}",
                f"samtools index {sorted_bam}"
            ]
            job_id = self.submit_slurm_job(cmds, f"samtools_{base}", dependencies=dependencies)
            jobs.append(job_id)

        return jobs

    def run_genome_coverage(self, dependencies: List[str]) -> List[str]:
        """Calcul de la couverture génomique avec bedtools genomecov"""
        jobs = []
        input_dir = Path(self.directories["sorted_bam"])
        output_dir = Path(self.directories["genomecov"])
        output_dir.mkdir(parents=True, exist_ok=True)

        for bam in input_dir.glob('*_sorted.bam'):
            bed = output_dir / f"{bam.stem}_genomecov.bed"
            if bed.exists():
                logging.info(f"Couverture {bed} existe déjà → skip")
                continue
            cmd = f"bedtools genomecov -ibam {bam} -bga > {bed}"
            job_id = self.submit_slurm_job([cmd], f"gencov_{bam.stem}", dependencies=dependencies)
            jobs.append(job_id)

        return jobs

    def extract_low_coverage(self, dependencies: List[str]) -> List[str]:
        """Extraction des régions à faible couverture (<30x)"""
        jobs = []
        input_dir = Path(self.directories["genomecov"])
        output_dir = Path(self.directories["extractline"])
        output_dir.mkdir(parents=True, exist_ok=True)

        for bed in input_dir.glob('*_genomecov.bed'):
            base = bed.stem.replace('_genomecov', '')
            out_bed = output_dir / f"{base}_low.bed"
            if out_bed.exists():
                logging.info(f"Fichier low coverage {out_bed} existe déjà → skip")
                continue
            cmd = f"awk '$4 <= 30 {{print}}' {bed} > {out_bed}"
            job_id = self.submit_slurm_job([cmd], f"extract_{base}", dependencies=dependencies)
            jobs.append(job_id)

        return jobs


    def call_variants(self, dependencies: List[str]) -> List[str]:
        """Appel des variants avec samtools mpileup + bcftools"""
        jobs = []
        input_dir = Path(self.directories["sorted_bam"])
        output_dir = Path(self.directories["variants"])
        output_dir.mkdir(parents=True, exist_ok=True)

        for bam in input_dir.glob('*_sorted.bam'):
            base = bam.stem.replace('_sorted', '')
            ref_key = 'rsv_a' if 'refA' in base else 'rsv_b'
            vcf = output_dir / f"{base}.vcf"

            if vcf.exists():
                logging.info(f"Fichier VCF {vcf} existe déjà → skip")
                continue

            # ✅ Nouvelle commande avec options valides
            cmd = (
                f"samtools mpileup -uf {self.references[ref_key]} {bam} "
                f"| bcftools call -mv -Ov -o {vcf}"
            )

            job_id = self.submit_slurm_job([cmd], f"mpileup_{base}", dependencies=dependencies)
            jobs.append(job_id)

        return jobs

    def compress_and_index_vcf(self, dependencies: List[str]) -> List[str]:
        """Compression et indexation des fichiers VCF"""
        jobs = []
        input_dir = Path(self.directories["variants"])

        for vcf in input_dir.glob('*.vcf'):
            vcf_gz = vcf.with_suffix('.vcf.gz')
            tbi_file = vcf.with_suffix('.vcf.gz.tbi')

            if vcf_gz.exists() and tbi_file.exists():
                logging.info(f"Fichiers {vcf_gz} et {tbi_file} existent déjà → skip")
                continue

            cmds = [
                f"bgzip -f {vcf}",
                f"bcftools index {vcf}.gz"
            ]
            job_id = self.submit_slurm_job(cmds, f"compress_{vcf.stem}", dependencies=dependencies)
            jobs.append(job_id)

        return jobs

    def generate_consensus(self, dependencies: List[str]) -> List[str]:
        """Génération du consensus à partir des variants et BED faible couverture"""
        jobs = []
        input_dir = Path(self.directories["variants"])
        extract_dir = Path(self.directories["extractline"])
        consensus_dir = Path(self.directories["consensus"])
        consensus_dir.mkdir(parents=True, exist_ok=True)

        for vcf in input_dir.glob('*.vcf.gz'):
            folder_name = vcf.parent.name
            base_part = vcf.stem
            base = base_part.split("_ref")[0]
            ref_type = "refA" if "_refA" in base_part else "refB"
            ref_key = 'rsv_a' if ref_type == 'refA' else 'rsv_b'

            bed_file = extract_dir / f"{base}_{ref_type}_sorted_low.bed"
            print(bed_file)
            consensus_file = consensus_dir / f"{base}_consensus_{ref_type}.fasta"

            if consensus_file.exists():
                logging.info(f"Consensus {consensus_file} existe déjà → skip")
                continue

            if not bed_file.exists():
                logging.warning(f"Fichier BED manquant pour {base} → consensus non lancé")
                continue

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
            logging.info("Début du pipeline Illumina RSV")
            trim_jobs = self.trim_adapters()
            bowtie_jobs = self.run_bowtie2(trim_jobs)
            sam_jobs = self.process_sam_files(bowtie_jobs)
            cov_jobs = self.run_genome_coverage(sam_jobs)
            extract_jobs = self.extract_low_coverage(cov_jobs)
            variant_jobs = self.call_variants(extract_jobs)
            compress_jobs = self.compress_and_index_vcf(variant_jobs)
            consensus_jobs = self.generate_consensus(compress_jobs)
            logging.info(f"Pipeline terminé. Consensus généré dans {self.directories['consensus']}")
        except Exception as e:
            logging.critical(f"Échec du pipeline: {str(e)}")
            raise


if __name__ == "__main__":
    pipeline = SlurmPipelineRSV()
    pipeline.execute_pipeline()