# Author: Dhatri Badri

# def samples_that_passed_coverage(wildcards=None, return_samples_only=False):
#     # Get the output file from the checkpoint
#     summary_csv = checkpoints.summarize_coverage.get(prefix=PREFIX).output[0]
#     # Read the coverage summary
#     df = pd.read_csv(summary_csv)
#     # Filter samples with coverage > 20
#     passed_samples = df[df['Coverage'] > 20]['Sample'].tolist()
#     if return_samples_only:
#         return passed_samples
#     return expand("results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html", prefix=PREFIX, sample=passed_samples) + expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html", prefix=PREFIX, sample=passed_samples) + expand("results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz", prefix=PREFIX, sample=passed_samples) 


checkpoint summarize_coverage:
    input:
        coverage_files = expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", sample=SAMPLE, prefix=PREFIX)
    output:
        "results/{prefix}/coverage_summary.csv"
    run:
        rows = []
        for cov_file in input.coverage_files:
            with open(cov_file) as f:
                data = json.load(f)
            sample = os.path.basename(cov_file).replace("_coverage.json", "")
            coverage = data.get('qc_stats', {}).get('coverage', 0)
            total_reads = data.get('qc_stats', {}).get('read_total', 0)
            total_bp = data.get('qc_stats', {}).get('total_bp', 0)
            mean_read_length = data.get('qc_stats', {}).get('read_mean', 0)
            rows.append({"Sample": sample, "Total_reads": total_reads, "Total_bp": total_bp, "MeanReadLength": mean_read_length, "Coverage": coverage})
        df = pd.DataFrame(rows)
        df.to_csv(output[0], index=False)