#!/bin/bash
#SBATCH --job-name=samtools_depth
#SBATCH --output=logs/samtools_depth_%j.out
#SBATCH --error=logs/samtools_depth_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --partition=fast

BAM_DIR="/shared/projects/virome_human_apes/Maud_biostats/pipeline_py_NanoporeAnalysis/data/output/sorted_bam"
DEPTH_DIR="$BAM_DIR/../depth"
mkdir -p "$DEPTH_DIR"

for bam in "$BAM_DIR"/*.bam; do
    [ -e "$bam" ] || continue
    base=$(basename "$bam" .bam)
    output_file="$DEPTH_DIR/${base}_depth.txt"
    
    echo "Processing $base..."
    samtools depth "$bam" > "$output_file"
done

echo "✅ Calcul de profondeur terminé."