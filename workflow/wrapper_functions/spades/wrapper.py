__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2025, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

from snakemake.shell import shell
import subprocess
import os
import pandas as pd
import json

out_dir = snakemake.params.get("out_dir", "")
threads = str(snakemake.params.get("threads", ""))
prefix_dir = snakemake.params.get("prefix_dir", "")
sample = snakemake.params.get("sample", "")

r1_files = snakemake.input.r1  # List from expand()
r2_files = snakemake.input.r2  # List from expand()

# print(f"r1_files looks like {r1_files}")

# # Ensure inputs are properly formatted for SPAdes
r1 = " ".join(r1_files) if isinstance(r1_files, list) else r1_files
r2 = " ".join(r2_files) if isinstance(r2_files, list) else r2_files
# print(f"r1 looks like {r1}")

spades_output_file = snakemake.output.spades_assembly  # Expected output file
#updated_samples_pass_assembly_successfully = snakemake.output.updated_samples_file_successsfully

def run_spades():
    """Runs SPAdes and captures stdout/stderr."""
    try:
        # assert r1, "Error: r1 input is missing"
        # assert r2, "Error: r2 input is missing"
        # assert out_dir, "Error: out_dir is missing"
        # assert threads, "Error: threads is missing"

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
        return result.stdout, None  # No errors

    except subprocess.CalledProcessError as e:
        return e.stdout + e.stderr, e.returncode  # Capture error output

def update_passed_csv(sample_name, outdir):
    """Update the CSV file for successfully assembled samples."""
    passed_csv = os.path.join(outdir, "sample_files", "samples_passed_assembly.csv") # passed_csv = os.path.join(outdir, "sample_files", "samples_passed_assembly_step.csv")
    
    # Read existing passed samples to avoid duplicates
    existing_passed_samples = set(pd.read_csv(passed_csv)["sample_id"])

    if sample_name not in existing_passed_samples:
        passed_df = pd.DataFrame([[sample_name, f"{sample_name}_R1.fastq.gz"]],
                                columns=["sample_id", "illumina_r1"])
        
        #header_needed = not os.path.exists(passed_csv)
        passed_df.to_csv(passed_csv, mode='a', index=False, header=False)

        # with open(updated_samples_pass_assembly_successfully, "a") as done_file:
        #     done_file.write(f"Sample: {sample_name} updated successfully.\n")

def update_failed_csv(sample_name, coverage_json, outdir):
    """Update the CSV file for failed assemblies."""
    failed_csv = os.path.join(outdir, "sample_files", "samples_failed_assembly_summary.csv")
    
    print(f"path to coverage file: {coverage_json}")

    try:
        with open(coverage_json) as f:
            data = json.load(f)
        coverage = data.get('qc_stats', {}).get('coverage', 0)
        total_reads = data.get('qc_stats', {}).get('read_total', 0)
        total_bp = data.get('qc_stats', {}).get('total_bp', 0)
        mean_read_length = data.get('qc_stats', {}).get('read_mean', 0)
    except (FileNotFoundError, json.JSONDecodeError):
        coverage, total_reads, total_bp, mean_read_length = 0, 0, 0, 0

    failed_df = pd.DataFrame(
        [[sample_name, total_reads, total_bp, mean_read_length, coverage] + ["NA"] * 15 + ["FAIL"] + ["NA"] * 5],
        columns=[
            "Sample", "Total_reads", "Total_bp", "MeanReadLength", "Coverage", "Scheme",
            "ST", "After_trim_per_base_sequence_content", "After_trim_overrepresented_sequences",
            "After_trim_%GC", "After_trim_Total Bases", "After_trim_Total Sequences",
            "After_trim_median_sequence_length", "After_trim_avg_sequence_length",
            "After_trim_total_deduplicated_percentage", "After_trim_Sequence length",
            "After_trim_adapter_content", "N50", "Total length", "Total # of contigs",
            "QC Check", "ANI", "Align_fraction_ref", "Align_fraction_query",
            "Ref_name", "Species"
        ]
    )
    
    header_needed = not os.path.exists(failed_csv)
    failed_df.to_csv(failed_csv, mode='a', index=False, header=header_needed)

def main():
    coverage_json = os.path.join(prefix_dir, "raw_coverage", sample, f"{sample}_coverage.json")  
    contigs_fasta = os.path.join(out_dir, "contigs.fasta") 
    
    # Run SPAdes and capture output
    spades_output, error_code = run_spades()

    # Check if the assembly failed due to OS return value 21
    if error_code == 21:
        update_failed_csv(sample, coverage_json, prefix_dir)
        print(f"Assembly failed for {sample} due to OS return value 21. Updated samples_failed_assembly_summary.csv.")

        # Create an empty contigs.fasta to prevent Snakemake from waiting on a missing assembly file
        with open(contigs_fasta, "w") as f:
            f.write(">failed_assembly\n")  # Add a dummy sequence header
        print(f"Created an empty contigs.fasta for failed sample {sample}.")
        print(f"Error code is: {error_code}")
    else:
        # If SPAdes ran successfully and the output file exists, mark the sample as passed
        if os.path.exists(spades_output_file):  
            update_passed_csv(sample, prefix_dir)
            print(f"Assembly successful for {sample}. Updated samples_passed_assembly.csv.") # print(f"Assembly successful for {sample}. Updated samples_passed_assembly_step.csv.")
        else:
            print(f"Error code is: {error_code}")
            update_failed_csv(sample, coverage_json, prefix_dir)
            with open(contigs_fasta, "w") as f:
                f.write(">failed_assembly\n")  # Add a dummy sequence header
            print(f"Created an empty contigs.fasta for failed sample {sample}.")
            print(f"Unexpected failure for {sample}, but OS return code was not 21. Logging output:")
            print(spades_output)

if __name__ == "__main__":
    main()
