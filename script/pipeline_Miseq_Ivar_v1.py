#!/usr/bin/env python3
import subprocess
import tempfile
import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, Any, List
from utils import setup_directories, setup_logging  # À adapter selon tes outils


class SlurmPipelineRSV:
    def __init__(self, config_path: str = "config/paths3.json"):
        self.config = self.load_config(config_path)
        self.references = self.config["references"]
        self.directories = self.config["data"]["output"]
        self.slurm_params = self.config["slurm"]
        self.modules = self.config["modules"]

        # Création des dossiers de sortie
        setup_directories(self.directories)

        # Configuration du logging
        setup_logging(self.config["logs"])

        logging.info("Initialisation du pipeline RSV Illumina avec Bowtie2 + ivar")

    def load_config(self, config_path: str) -> Dict[str, Any]:
        with open(config_path) as f:
            return json.load(f)

    def create_slurm_script(self, commands: List[str], job_name: str, dependencies: List[str] = None) -> str:
        """Crée un script Slurm temporaire"""

        # Générer les lignes pour charger les modules
        module_commands = [f"module load {mod}" for mod in self.modules.values()]
        module_lines = "\n".join(module_commands)

        # Ligne pour les dépendances (si présentes)
        dependency_line = ""
        if dependencies:
            dependency_line = f"#SBATCH --dependency=afterok:{':'.join(dependencies)}"

        # Préparer les commandes à exécuter
        command_lines = "\n".join(commands)

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
    {command_lines}
    """

        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(slurm_template)
            return f.name

    def submit_slurm_job(self, commands: List[str], step_name: str, dependencies: List[str] = None) -> str:
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
        jobs = []
        input_dir = Path(self.config["data"]["input"])
        output_dir = Path(self.directories["trimmed"])
        output_dir.mkdir(parents=True, exist_ok=True)

        for r1 in input_dir.glob('*_R1*.fastq.gz'):
            base = r1.name.replace('_R1', '').replace('.fastq.gz', '')
            r2 = input_dir / r1.name.replace('_R1', '_R2')
            if not r2.exists():
                logging.warning(f"Fichier R2 manquant pour {base} → skip")
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

    def generate_consensus(self, dependencies: List[str]) -> List[str]:
        """Génère le consensus directement avec ivar consensus depuis le BAM trié"""
        jobs = []
        input_dir = Path(self.directories["sorted_bam"])
        consensus_dir = Path(self.directories["consensus"])
        consensus_dir.mkdir(parents=True, exist_ok=True)

        for bam in input_dir.glob('*_sorted.bam'):
            base = bam.stem.replace('_sorted', '')  # ex: sample_refA
            ref_type = "refA" if "_refA" in base else "refB"
            ref_key = 'rsv_a' if ref_type == 'refA' else 'rsv_b'
            consensus_file = consensus_dir / f"{base}_ivar_consensus.fasta"

            if consensus_file.exists():
                logging.info(f"Consensus {consensus_file} existe déjà → skip")
                continue

            # Commande complète avec activation de l'environnement conda
            cmd = (
                f"source activate ivar_env && "
                f"samtools mpileup -aa -A -d 0 -Q 0 {bam} | "
                f"ivar consensus -p  {consensus_file}"
            )

            job_id = self.submit_slurm_job([cmd], f"ivar_{base}", dependencies=dependencies)
            jobs.append(job_id)

        return jobs
    def execute_pipeline(self):
        try:
            logging.info("Début du pipeline Illumina RSV")
            trim_jobs = self.trim_adapters()
            bowtie_jobs = self.run_bowtie2(trim_jobs)
            sam_jobs = self.process_sam_files(bowtie_jobs)
            consensus_jobs = self.generate_consensus(sam_jobs)
            logging.info(f"Pipeline terminé. Consensus généré dans {self.directories['consensus']}")
        except Exception as e:
            logging.critical(f"Échec du pipeline : {str(e)}")
            raise


if __name__ == "__main__":
    pipeline = SlurmPipelineRSV()
    pipeline.execute_pipeline()