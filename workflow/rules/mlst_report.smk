rule mlst_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        mlst_out = lambda wildcards: [
            f"results/{PREFIX}/mlst/{sample}/report.tsv"
            for sample in samples_that_passed_assembly(wildcards, return_samples_only=True)
        ] # expand("results/{prefix}/mlst/{sample}/report.tsv", prefix=PREFIX, sample=SAMPLE)
    output:
        mlst_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_MLST_results.csv",
    params:
        prefix = "{prefix}",
    shell:
        "echo \"Sample\tScheme\tST\" > {output.mlst_report} && cut -f1-3 {input.outdir}/mlst/*/report.tsv >> {output.mlst_report}"
