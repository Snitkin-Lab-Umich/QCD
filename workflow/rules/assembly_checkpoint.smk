# def samples_that_passed_assembly(wildcards=None, return_samples_only=False):
#     # Get the output file from the checkpoint
#     summary_csv = checkpoints.summarize_assembly.get(prefix=PREFIX).output[0]
#     # Read the coverage summary
#     df = pd.read_csv(summary_csv)
#     # Filter samples where Spades_success is True
#     passed_samples_assembly = df[df['Spades_success'] == True]['Sample'].tolist()
#     if return_samples_only:
#         return passed_samples_assembly
#     return expand("results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta", prefix=PREFIX, sample=passed_samples_assembly) + \
#             expand("results/{prefix}/prokka/{sample}/{sample}.gff", prefix=PREFIX, sample=passed_samples_assembly) + \
#             expand("results/{prefix}/quast/{sample}/report.tsv", prefix=PREFIX, sample=passed_samples_assembly) + \
#             expand("results/{prefix}/mlst/{sample}/report.tsv", prefix=PREFIX, sample=passed_samples_assembly) + \
#             expand("results/{prefix}/busco/{sample}/{sample}_busco_out.txt", prefix=PREFIX, sample=passed_samples_assembly) + \
#             expand("results/{prefix}/skani/{sample}/{sample}_skani_output.txt", prefix=PREFIX, sample=passed_samples_assembly) 

checkpoint summarize_assembly:
    input:
        spades_assemblies = lambda wildcards: [
            f"results/{PREFIX}/spades/{sample}/contigs.fasta"
            for sample in samples_that_passed_coverage(wildcards, return_samples_only=True)
        ]
    output:
        "results/{prefix}/assembly_summary.csv"
    run:
        rows = []
        for spades_file in input.spades_assemblies:
            sample = os.path.basename(os.path.dirname(spades_file))
            spades_success = os.path.getsize(spades_file) > 0  # not empty
            rows.append({"Sample": sample, "Spades_success": spades_success})
        df = pd.DataFrame(rows)
        df.to_csv(output[0], index=False)