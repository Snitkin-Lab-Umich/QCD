rule multiqc:
    input:
        inputdir = lambda wildcards: expand(f"results/{wildcards.prefix}"),
        coverage = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Final_Coverage.txt"),
        #kraken = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Kraken_report_final.csv"),
        mlst = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_MLST_results.csv"),
    output:
        multiqc_fastqc_report = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report.html",
        multiqc_fastqc = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report_data/multiqc_fastqc.txt",
        multiqc_general_stats = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report_data/multiqc_general_stats.txt",
        #fastqc_report = f"results//{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report_data/multiqc_fastqc.txt"
    params:
        outdir = "results/{prefix}/{prefix}_Report",
        prefix = "{prefix}",
    #conda:
    #    "envs/multiqc.yaml"
    singularity:
        "docker://staphb/multiqc:1.19"
    shell:
        "multiqc -f --export --outdir {params.outdir}/multiqc -n {params.prefix}_QC_report -i {params.prefix}_QC_report {input.inputdir}/quality_aftertrim/*/*_Forward {input.inputdir}/prokka/* {input.inputdir}/quast/*"
