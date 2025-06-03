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
        setup_logging(self.config["logs"])
        logging.info("Initialisation du pipeline RSV avec Slurm")

    def load_config(self, config_path: str) -> Dict[str, Any]:
        """Charge la configuration depuis le fichier JSON"""
        with open(config_path) as f:
            return json.load(f)

    def create_slurm_script(self, commands: List[str], job_name: str, dependencies: List[str] = None) -> str:
        """Crée un script Slurm temporaire"""
        slurm_template = f"""#!/bin/bash
#SBATCH --account={self.slurm_params["account"]}
#SBATCH --job-name={job_name}
#SBATCH --output={self.config["logs"]}/{job_name}_%j.out
#SBATCH --error={self.config["logs"]}/{job_name}_%j.err
#SBATCH --time={self.slurm_params["time"]}
#SBATCH --mem={self.slurm_params["mem"]}
#SBATCH --cpus-per-task={self.slurm_params["cpus"]}
#SBATCH --partition={self.slurm_params["partition"]}

{dependencies and f"#SBATCH --dependency=afterok:{':'.join(dependencies)}" or ''}

# Chargement des modules
module purge
{"".join(f"module load {mod}\n" for mod in self.modules.values())}

# Commandes
{"".join(f"{cmd}\n" for cmd in commands)}
"""

        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(slurm_template)
            return f.name
#-------------------------------------------------------------------------------------------------------------------------
   
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
#-------------------------------------------------------------------------------------------------------------------------
    def concatenate_fastq(self) -> List[str]:
        jobs = []
        input_dir = Path(self.config["data"]["input"])
        output_dir = Path(self.directories["concatenated"])
        
        # Créer le dossier de sortie si inexistant
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
        
        # Créer le dossier de sortie si inexistant
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Générer la liste des fichiers attendus basée sur l'input initial
        input_parent = Path(self.config["data"]["input"])
        expected_files = [input_dir / f"{d.name}.fastq.gz" 
                        for d in input_parent.iterdir() 
                        if d.is_dir()]
        
        for fq in expected_files:
            base = fq.stem  # Retire l'extension .fastq.gz
            output_file = output_dir / f"{base}_trimmed.fastq.gz"
            
            if output_file.exists():
                logging.info(f"Fichier trimmed {output_file} existe déjà - skip")
                continue
                
            cmd = f"cutadapt -u 20 -u -20 -o {output_file} {fq}"
            job_id = self.submit_slurm_job(
                [cmd], 
                f"cutadapt_{base}", 
                dependencies=concat_job_ids  # Dépendance sur tous les jobs de concat
            )
            jobs.append(job_id)
        
        return jobs
        
    def run_nanoplot(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_parent = Path(self.config["data"]["input"])
        output_parent = Path(self.directories["nanoplot"])
        trimmed_dir = Path(self.directories["trimmed"])
        
        output_parent.mkdir(parents=True, exist_ok=True)
        
        # Générer la liste des fichiers attendus basée sur la structure d'entrée
        expected_files = [trimmed_dir / f"{d.name}_trimmed.fastq.gz" 
                        for d in input_parent.iterdir() 
                        if d.is_dir()]
        
        for fq in expected_files:
            base = fq.stem.replace('_trimmed', '')
            output_dir = output_parent / f"{base}_quality"
            
            if (output_dir / "NanoPlot-report.html").exists():
                logging.info(f"NanoPlot {output_dir} existe déjà - skip")
                continue
                
            cmd = f"NanoPlot --fastq {fq} --outdir {output_dir}"
            job_id = self.submit_slurm_job([cmd], f"nanoplot_{base}", dependencies)
            jobs.append(job_id)
        
        return jobs

    def run_minimap2(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_parent = Path(self.config["data"]["input"])
        output_dir = Path(self.directories["aligned"])
        trimmed_dir = Path(self.directories["trimmed"])
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Générer la liste des fichiers attendus
        expected_files = [trimmed_dir / f"{d.name}_trimmed.fastq.gz" 
                        for d in input_parent.iterdir() 
                        if d.is_dir()]
        
        for fq in expected_files:
            base = fq.stem.replace('_trimmed', '')
            
            # RSV-A
            sam_a = output_dir / f"{base}_refA.sam"
            cmd_a = f"minimap2 -ax map-ont {self.references['rsv_a']} {fq} > {sam_a}"
            job_id = self.submit_slurm_job([cmd_a], f"minimapA_{base}", dependencies)
            jobs.append(job_id)
            
            # RSV-B
            sam_b = output_dir / f"{base}_refB.sam"
            cmd_b = f"minimap2 -ax map-ont {self.references['rsv_b']} {fq} > {sam_b}"
            job_id = self.submit_slurm_job([cmd_b], f"minimapB_{base}", dependencies)
            jobs.append(job_id)
        
        return jobs

    def process_sam_files(self, dependencies: List[str]) -> List[str]:
        jobs = []
        input_parent = Path(self.config["data"]["input"])
        output_dir = Path(self.directories["sorted_bam"])
        aligned_dir = Path(self.directories["aligned"])
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Générer la liste des SAM attendus
        expected_sams = []
        for d in input_parent.iterdir():
            if d.is_dir():
                base = d.name
                expected_sams.append(aligned_dir / f"{base}_refA.sam")
                expected_sams.append(aligned_dir / f"{base}_refB.sam")
        
        for sam in expected_sams:
            base = sam.stem
            bam = output_dir / f"{base}_sorted.bam"
            
            if bam.exists():
                logging.info(f"BAM {bam} existe déjà - skip")
                continue
                
            cmds = [
                f"samtools view -b {sam} | samtools sort -o {bam}",
                f"samtools index {bam}"
            ]
            job_id = self.submit_slurm_job(cmds, f"samtools_{base}", dependencies)
            jobs.append(job_id)
        
        return jobs

    def execute_pipeline(self):
        try:
            logging.info("Début du pipeline")
            
            # Chaînage séquentiel
            concat_jobs = self.concatenate_fastq()
            trim_jobs = self.trim_adapters(concat_jobs)
            nanoplot_jobs = self.run_nanoplot(trim_jobs)
            minimap_jobs = self.run_minimap2(trim_jobs)
            samtools_jobs = self.process_sam_files(minimap_jobs)
            
            logging.info(f"Jobs finals soumis: {samtools_jobs}")
        except Exception as e:
            logging.critical(f"Échec du pipeline: {str(e)}")
            raise

if __name__ == "__main__":
    pipeline = SlurmPipelineRSV()
    pipeline.execute_pipeline()
    
    def run_nanoplot(self, dependencies: List[str]) -> List[str]:
        
        jobs = []
        input_dir = Path(self.directories["trimmed"])
        
        for fq in input_dir.glob('*_trimmed.fastq.gz'):
            base = fq.name.replace('_trimmed.fastq.gz', '')
            output_dir = Path(self.directories["nanoplot"]) / f"{base}_quality"
            cmd = f"NanoPlot --fastq {fq} --outdir {output_dir}"
            job_id = self.submit_slurm_job([cmd], f"nanoplot_{base}", dependencies)
            jobs.append(job_id)
        
        return jobs

    def run_minimap2(self, dependencies: List[str]) -> List[str]:
        
        jobs = []
        input_dir = Path(self.directories["trimmed"])
        output_dir = Path(self.directories["aligned"])
        
        for fq in input_dir.glob('*_trimmed.fastq.gz'):
            # Alignement RSV-A
            base = fq.name.replace('_trimmed.fastq.gz', '')
            sam_a = output_dir / f"{base}_refA.sam"
            cmd_a = f"minimap2 -a {self.references['rsv_a']} {fq} > {sam_a}"
            job_id_a = self.submit_slurm_job([cmd_a], f"minimapA_{base}", dependencies)
            jobs.append(job_id_a)
            
            # Alignement RSV-B
            sam_b = output_dir / f"{base}_refB.sam"
            cmd_b = f"minimap2 -a {self.references['rsv_b']} {fq} > {sam_b}"
            job_id_b = self.submit_slurm_job([cmd_b], f"minimapB_{base}", dependencies)
            jobs.append(job_id_b)
        
        return jobs

    def process_sam_files(self, dependencies: List[str]) -> List[str]:
        
        jobs = []
        input_dir = Path(self.directories["aligned"])
        output_dir = Path(self.directories["sorted_bam"])
        
        for sam in input_dir.glob('*_refA.sam'):
            # Alignement RSV-A
            base = sam.name.replace('_refA.sam', '')
            bam_a = output_dir / f"{base}_refA_sorted.bam"
            cmd_a = [
                f"samtools view -b {sam} | samtools sort -o {bam_a}",
                f"samtools index {bam_a}"
            ]
            job_id_a = self.submit_slurm_job([cmd_a], f"samtoolsA_{base}", dependencies)
            jobs.append(job_id_a)
            
            # Alignement RSV-B
            sam_b = input_dir / f"{base}_refB.sam"
            bam_b = output_dir / f"{base}_refB_sorted.bam"
            cmd_b = [
                f"samtools view -b {sam_b} | samtools sort -o {bam_b}",
                f"samtools index {bam_b}"
            ]
            job_id_b = self.submit_slurm_job([cmd_b], f"samtoolsB_{base}", dependencies)
            jobs.append(job_id_b)
        
        return jobs

    def run_genome_coverage(self, dependencies: List[str]) -> List[str]:
        
        jobs = []
        input_dir = Path(self.directories["sorted_bam"])
        output_dir = Path(self.directories["genomecov"])
        
        for bam in input_dir.glob('*_refA_sorted.bam'):
            base = bam.name.replace('_refA_sorted.bam', '')
            # Couverture RSV-A
            bed_a = output_dir / f"{base}_refA_genomecov.bed"
            cmd_a = f"bedtools genomecov -ibam {bam} -bga > {bed_a}"
            job_id_a = self.submit_slurm_job([cmd_a], f"gencovA_{base}", dependencies)
            jobs.append(job_id_a)
            
            # Couverture RSV-B
            bam_b = input_dir / f"{base}_refB_sorted.bam"
            bed_b = output_dir / f"{base}_refB_genomecov.bed"
            cmd_b = f"bedtools genomecov -ibam {bam_b} -bga > {bed_b}"
            job_id_b = self.submit_slurm_job([cmd_b], f"gencovB_{base}", dependencies)
            jobs.append(job_id_b)
        
        return jobs

    def extract_low_coverage(self, dependencies: List[str]) -> List[str]:
        
        jobs = []
        input_dir = Path(self.directories["genomecov"])
        output_dir = Path(self.directories["extractline"])
        
        for bed in input_dir.glob('*_refA_genomecov.bed'):
            base = bed.name.replace('_refA_genomecov.bed', '')
            # RSV-A
            cmd_a = f"awk '$4 <= 30 {{print}}' {bed} > {output_dir}/{base}_refA_extractline.bed"
            job_id_a = self.submit_slurm_job([cmd_a], f"extractA_{base}", dependencies)
            jobs.append(job_id_a)
            
            # RSV-B
            bed_b = input_dir / f"{base}_refB_genomecov.bed"
            cmd_b = f"awk '$4 <= 30 {{print}}' {bed_b} > {output_dir}/{base}_refB_extractline.bed"
            job_id_b = self.submit_slurm_job([cmd_b], f"extractB_{base}", dependencies)
            jobs.append(job_id_b)
        
        return jobs

    def call_variants(self, dependencies: List[str]) -> List[str]:
        
        jobs = []
        input_dir = Path(self.directories["sorted_bam"])
        output_dir = Path(self.directories["medaka_variants"])
        
        for bam in input_dir.glob('*_refA_sorted.bam'):
            base = bam.name.replace('_refA_sorted.bam', '')
            
            # RSV-A
            output_a = output_dir / f"{base}_medaka_variant_refA"
            cmd_a = f"medaka_variant -i {bam} -r {self.references['rsv_a']} -o {output_a}"
            job_id_a = self.submit_slurm_job([cmd_a], f"medakaA_{base}", dependencies)
            jobs.append(job_id_a)
            
            # RSV-B
            bam_b = input_dir / f"{base}_refB_sorted.bam"
            output_b = output_dir / f"{base}_medaka_variant_refB"
            cmd_b = f"medaka_variant -i {bam_b} -r {self.references['rsv_b']} -o {output_b}"
            job_id_b = self.submit_slurm_job([cmd_b], f"medakaB_{base}", dependencies)
            jobs.append(job_id_b)
        
        return jobs

    def compress_and_index_vcf(self, dependencies: List[str]) -> List[str]:
       
        jobs = []
        input_dir = Path(self.directories["medaka_variants"])
        
        for vcf in input_dir.rglob('medaka.sorted.vcf'):
            cmds = [
                f"bgzip {vcf}",
                f"bcftools index {vcf}.gz"
            ]
            job_id = self.submit_slurm_job(cmds, f"compress_{vcf.parent.name}", dependencies)
            jobs.append(job_id)
        
        return jobs

    def generate_consensus(self, dependencies: List[str]) -> List[str]:
        
        jobs = []
        input_dir = Path(self.directories["medaka_variants"])
        extract_dir = Path(self.directories["extractline"])
        consensus_dir = Path(self.directories["consensus"])
        
        for ref_type in ['refA', 'refB']:
            for vcf in input_dir.glob(f'*_medaka_variant_{ref_type}/medaka.sorted.vcf.gz'):
                base = vcf.parent.name.replace(f'_medaka_variant_{ref_type}', '')
                ref_key = 'rsv_a' if ref_type == 'refA' else 'rsv_b'
                bed_file = extract_dir / f"{base}_{ref_type}_extractline.bed"
                consensus_file = consensus_dir / f"{base}_consensus_{ref_type}.fasta"
                
                cmd = (
                    f"bcftools consensus -m {bed_file} "
                    f"-f {self.references[ref_key]} -o {consensus_file} {vcf}"
                )
                job_id = self.submit_slurm_job([cmd], f"consensus_{ref_type}_{base}", dependencies)
                jobs.append(job_id)
        
        return jobs

            # nanoplot_jobs = self.run_nanoplot(trim_jobs)
            # minimap_jobs = self.run_minimap2(nanoplot_jobs)
            # sam_jobs = self.process_sam_files(minimap_jobs)
            # cov_jobs = self.run_genome_coverage(sam_jobs)
            # extract_jobs = self.extract_low_coverage(cov_jobs)
            # medaka_jobs = self.call_variants(extract_jobs)
            # compress_jobs = self.compress_and_index_vcf(medaka_jobs)
            # consensus_jobs = self.generate_consensus(compress_jobs)