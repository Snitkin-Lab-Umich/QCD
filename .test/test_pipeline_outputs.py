import pandas as pd
import os

# Helper function to compare CSVs ignoring row order and whitespace
def compare_csv_files(file1, file2):
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    # Compare shapes
    if df1.shape != df2.shape:
        return False
    # Sort by all columns (assuming no index column)
    df1_sorted = df1.sort_values(by=df1.columns.tolist()).reset_index(drop=True)
    df2_sorted = df2.sort_values(by=df2.columns.tolist()).reset_index(drop=True)
    # Compare all values
    return df1_sorted.equals(df2_sorted)

def test_skani_report_output():
    prefix = "2025-05-29_Project_MRSA_QCD"  # Prefix used in config.yaml
    current_file = f"results/{prefix}/{prefix}_Report/data/{prefix}_Skani_report_final.csv"
    reference_file = f".test/reference_data/{prefix}_Skani_report_final.csv"
    assert os.path.exists(current_file), f"{current_file} does not exist."
    assert os.path.exists(reference_file), f"{reference_file} does not exist."
    assert compare_csv_files(current_file, reference_file), "Skani report output differs from reference!"

def test_qc_summary_output():
    prefix = "2025-05-29_Project_MRSA_QCD"
    current_file = f"results/{prefix}/{prefix}_Report/data/{prefix}_QC_summary.csv"
    reference_file = f".test/reference_data/{prefix}_QC_summary.csv"
    assert os.path.exists(current_file), f"{current_file} does not exist."
    assert os.path.exists(reference_file), f"{reference_file} does not exist."
    assert compare_csv_files(current_file, reference_file), "QC summary output differs from reference!"

def test_mlst_report_output():
    prefix = "2025-05-29_Project_MRSA_QCD"
    current_file = f"results/{prefix}/{prefix}_Report/data/{prefix}_MLST_results.csv"
    reference_file = f".test/reference_data/{prefix}_MLST_results.csv"
    assert os.path.exists(current_file), f"{current_file} does not exist."
    assert os.path.exists(reference_file), f"{reference_file} does not exist."
    assert compare_csv_files(current_file, reference_file), "MLST report output differs from reference!"

# def test_quast_report_output():
#     prefix = "2025-05-29_Project_MRSA_QCD"
#     samples = ["MRSA_jail_100"]  # Sample that passed
#     for sample in samples:
#         current_file = f"results/{prefix}/quast/{sample}/report.tsv"
#         reference_file = f".test/reference_data/{sample}_report.tsv"
#         assert os.path.exists(current_file), f"{current_file} does not exist."
#         assert os.path.exists(reference_file), f"{reference_file} does not exist."
#         df_current = pd.read_csv(current_file, sep="\t")
#         df_ref = pd.read_csv(reference_file, sep="\t")
#         assert df_current.equals(df_ref), f"Quast report for {sample} differs from reference!"
