__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2025, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

from snakemake.shell import shell
import subprocess
import os
import pandas as pd
import json
import csv

out_dir = snakemake.params.get("out_dir", "")
threads = str(snakemake.params.get("threads", ""))
prefix_dir = snakemake.params.get("prefix_dir", "")
sample = snakemake.params.get("sample", "")

r1_files = snakemake.input.r1
r2_files = snakemake.input.r2

r1 = " ".join(r1_files) if isinstance(r1_files, list) else r1_files
r2 = " ".join(r2_files) if isinstance(r2_files, list) else r2_files

spades_output_file = snakemake.output.spades_assembly
summary_csv = os.path.join(out_dir, "assembly_summary.csv")

def run_spades():
    try:
        result = subprocess.run(
            [
                "spades.py", "--isolate",
                "--pe1-1", r1,
                "--pe1-2", r2,
                "-o", out_dir,
                "--threads", str(threads)
            ],
            capture_output=True, text=True, check=True
        )
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        return False, e.stdout + e.stderr

# def get_coverage_stats(coverage_json):
#     try:
#         with open(coverage_json) as f:
#             data = json.load(f)
#         coverage = data.get('qc_stats', {}).get('coverage', 0)
#         total_reads = data.get('qc_stats', {}).get('read_total', 0)
#         total_bp = data.get('qc_stats', {}).get('total_bp', 0)
#         mean_read_length = data.get('qc_stats', {}).get('read_mean', 0)
#     except Exception:
#         coverage, total_reads, total_bp, mean_read_length = 0, 0, 0, 0
#     return coverage, total_reads, total_bp, mean_read_length

# def write_assembly_summary(summary_csv, success, coverage, total_reads, total_bp, mean_read_length):
#     with open(summary_csv, "w", newline="") as csvfile:
#         writer = csv.writer(csvfile)
#         writer.writerow(["sample", "spades_success", "coverage", "total_reads", "total_bp", "mean_read_length"])
#         writer.writerow([
#             sample,
#             "yes" if success else "no",
#             coverage,
#             total_reads,
#             total_bp,
#             mean_read_length
#         ])

def main():
    coverage_json = os.path.join(prefix_dir, "raw_coverage", sample, f"{sample}_coverage.json")
    contigs_fasta = os.path.join(out_dir, "contigs.fasta")

    spades_success, spades_output = run_spades()
    # coverage, total_reads, total_bp, mean_read_length = get_coverage_stats(coverage_json)
    # write_assembly_summary(summary_csv, spades_success, coverage, total_reads, total_bp, mean_read_length)

    if not spades_success or not os.path.exists(spades_output_file):
        # Create empty contigs.fasta if failed
        with open(contigs_fasta, "w") as f:
            pass  # truly empty file
        print(f"SPAdes failed for {sample}. Created empty contigs.fasta.")
    else:
        print(f"SPAdes succeeded for {sample}.")

if __name__ == "__main__":
    main()