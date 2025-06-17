# Author: Dhatri Badri

rule quality_raw:
    input:
        # r1_trimmed = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"), # to ensure this runs aftrer DAG evaluates
        r1 = lambda wildcards: os.path.join(SHORT_READS_DIR, f"{wildcards.sample}_R1.fastq.gz"),
        r2 = lambda wildcards: os.path.join(SHORT_READS_DIR, f"{wildcards.sample}_R2.fastq.gz"),
    output:
        raw_fastqc_report_fwd = "results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html",
        raw_fastqc_report_rev = "results/{prefix}/quality_raw/{sample}/{sample}_Reverse/{sample}_R2_fastqc.html",
    log:
        "logs/{prefix}/quality_raw/{sample}/{sample}.log"
    params:
        outdir = "results/{prefix}/quality_raw/{sample}/{sample}"
    singularity:
        "docker://staphb/fastqc:0.12.1"
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        """
