#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Génère un résumé statistique à partir des fichiers .depth.txt")
    parser.add_argument("--input_dir", "-i", type=str, required=True, help="Dossier contenant les fichiers .depth.txt")
    parser.add_argument("--output", "-o", type=str, default="depth_summary.csv", help="Fichier CSV de sortie")
    parser.add_argument("--min_depth", "-m", type=int, default=30, help="Seuil de profondeur considéré comme bon")
    return parser.parse_args()

def summarize_depth_file(file_path, min_depth):
    depths = []
    covered_positions = 0
    total_positions = 0

    with open(file_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            depth = int(parts[2])
            depths.append(depth)
            if depth >= min_depth:
                covered_positions += 1
            total_positions += 1

    if not depths:
        return None

    mean_depth = np.mean(depths)
    median_depth = np.median(depths)
    coverage_percent = (covered_positions / total_positions) * 100 if total_positions > 0 else 0

    return {
        "sample": os.path.basename(file_path).replace("_depth.txt", ""),
        "total_positions": total_positions,
        "mean_depth": round(mean_depth, 2),
        "median_depth": round(median_depth, 2),
        "coverage_geq_min_depth": round(coverage_percent, 2),
        "positions_geq_min_depth": covered_positions
    }

def main():
    args = parse_args()

    depth_files = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) if f.endswith("_depth.txt")]

    results = []
    for file in depth_files:
        logging.info(f"Traitement de {file}")
        summary = summarize_depth_file(file, args.min_depth)
        if summary:
            results.append(summary)

    df = pd.DataFrame(results)
    df.to_csv(args.output, index=False)
    print(f"✅ Résumé écrit dans {args.output}")

if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO)
    main()