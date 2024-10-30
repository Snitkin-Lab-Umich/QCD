import os
import pytest
import subprocess

# Define the sample list
SAMPLE_LIST = [
    "INT_CRE_1196",
]

# Base directory for your results
RESULTS_DIR = "results/2024-10-22_11414-ES_run5_Plate14_QCD"


def run_snakemake():
    """Run the Snakemake pipeline."""
    command = ["snakemake", "--cores", "3", "--snakefile", "QCD.smk", "--dryrun"]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert result.returncode == 0, f"Snakemake failed with error: {result.stderr.decode()}"

@pytest.mark.parametrize("sample", SAMPLE_LIST)
def test_results_directory_structure():
    """Test if the output directory structure is created as expected."""
    # Base directory for your results
    base_dir = RESULTS_DIR

    # Extract the folder name
    qcd_folder_name = os.path.basename(RESULTS_DIR)

    # List of expected directories
    expected_dirs = [
        "downsample",
        "mlst",
        "prokka",
        "quality_aftertrim",
        "quality_raw",
        "quast",
        "raw_coverage",
        "skani",
        "spades",
        "trimmomatic",
    ]
    
    # Check if expected directories exist
    for dir_name in expected_dirs:
        assert os.path.isdir(os.path.join(base_dir, dir_name)), f"Expected directory {dir_name} does not exist."

    for sample in SAMPLE_LIST:
        # Check specific files in the 'prokka' directory
        prokka_files = [
            f"{sample}.err",
            f"{sample}.ffn",
            f"{sample}.fsa",
            f"{sample}.gff",
            f"{sample}.sqn",
            f"{sample}.tsv",
            f"{sample}.faa",
            f"{sample}.fna",
            f"{sample}.gbk",
            f"{sample}.log",
            f"{sample}.tbl",
            f"{sample}.txt"
        ]
        assert os.path.isdir(os.path.join(base_dir, "prokka", sample)), f"Expected directory for prokka does not exist for sample {sample}."
        for file_name in prokka_files:
            assert os.path.isfile(os.path.join(base_dir, "prokka", sample, file_name)), f"Expected file {file_name} does not exist in prokka directory for sample {sample}."

        # Check specific files in the downsample directory
        downsample_files = [
            f"{sample}_R1_trim_paired.fastq.gz",
            f"{sample}_R2_trim_paired.fastq.gz"
        ]
        assert os.path.isdir(os.path.join(base_dir, "downsample", sample)), f"Expected directory for downsample does not exist for sample {sample}."
        for file_name in downsample_files:
            assert os.path.isfile(os.path.join(base_dir, "downsample", sample, file_name)), f"Expected file {file_name} does not exist in downsample directory for sample {sample}."

        # Check specific files in the mlst directory
        assert os.path.isdir(os.path.join(base_dir, "mlst", sample)), f"Expected directory for mlst does not exist for sample {sample}."
        assert os.path.isfile(os.path.join(base_dir, "mlst", sample, "report.tsv")), f"Expected report.tsv does not exist in mlst directory for sample {sample}."

        # Check specific files in the quality_aftertrim directory
        assert os.path.isdir(os.path.join(base_dir, "quality_aftertrim", sample)), f"Expected directory for quality_aftertrim does not exist for sample {sample}."
        assert os.path.isfile(os.path.join(base_dir, "quality_aftertrim", sample, f"{sample}_Forward/{sample}_R1_trim_paired_fastqc.html")), f"Expected FastQC HTML file does not exist in quality_aftertrim directory for sample {sample}."

        # Check specific files in the quality_raw directory
        assert os.path.isdir(os.path.join(base_dir, "quality_raw", sample)), f"Expected directory for quality_raw does not exist for sample {sample}."
        assert os.path.isfile(os.path.join(base_dir, "quality_raw", sample, f"{sample}_Forward/{sample}_R1_fastqc.html")), f"Expected FastQC HTML file does not exist in quality_raw directory for sample {sample}."

        # Check specific files in the quast directory
        quast_files = [
            "report.html",
            "report.tsv",
            "report.txt"
        ]
        assert os.path.isdir(os.path.join(base_dir, "quast", sample)), f"Expected directory for quast does not exist for sample {sample}."
        for file_name in quast_files:
            assert os.path.isfile(os.path.join(base_dir, "quast", sample, file_name)), f"Expected file {file_name} does not exist in quast directory for sample {sample}."

        # Check specific files in the raw_coverage directory
        assert os.path.isdir(os.path.join(base_dir, "raw_coverage", sample)), f"Expected directory for raw_coverage does not exist for sample {sample}."
        assert os.path.isfile(os.path.join(base_dir, "raw_coverage", sample, f"{sample}_coverage.json")), f"Expected coverage.json does not exist in raw_coverage directory for sample {sample}."

        # Check specific files in the skani directory
        assert os.path.isdir(os.path.join(base_dir, "skani", sample)), f"Expected directory for skani does not exist for sample {sample}."
        assert os.path.isfile(os.path.join(base_dir, "skani", sample, f"{sample}_skani_output.txt")), f"Expected skani output file does not exist for sample {sample}."

        # Check specific files in the spades directory
        spades_files = [
            "contigs.fasta",
            f"{sample}_contigs_l1000.fasta"
        ]
        assert os.path.isdir(os.path.join(base_dir, "spades", sample)), f"Expected directory for spades does not exist for sample {sample}."
        for file_name in spades_files:
            assert os.path.isfile(os.path.join(base_dir, "spades", sample, file_name)), f"Expected file {file_name} does not exist in spades directory for sample {sample}."

        # Check specific files in the 'trimmomatic' directory
        trimmomatic_files = [
            f"{sample}_R1_trim_paired.fastq.gz",
            f"{sample}_R1_trim_unpaired.fastq.gz",
            f"{sample}_R2_trim_paired.fastq.gz",
            f"{sample}_R2_trim_unpaired.fastq.gz"
        ]
        assert os.path.isdir(os.path.join(base_dir, "trimmomatic", sample)), f"Expected directory for trimmomatic does not exist for sample {sample}."
        for file_name in trimmomatic_files:
            assert os.path.isfile(os.path.join(base_dir, "trimmomatic", sample, file_name)), f"Expected file {file_name} does not exist in trimmomatic directory for sample {sample}."

        # Check the report directory and its contents
        report_dir = os.path.join(base_dir, f"{qcd_folder_name}_Report", "data")
        assert os.path.isdir(report_dir), f"Expected report directory does not exist for {qcd_folder_name}."

        # List of expected files in the report data directory
        report_data_files = [
            f"{sample}_Final_Coverage.txt",
            f"{sample}_QC_summary.csv",
            f"{sample}_MLST_results.csv",
            f"{sample}_Skani_report_final.csv"
        ]
        for file_name in report_data_files:
            assert os.path.isfile(os.path.join(report_dir, file_name)), f"Expected file {file_name} does not exist in report data directory."

if __name__ == "__main__":
    run_snakemake()  # Run the pipeline to prepare for testing
    pytest.main()
