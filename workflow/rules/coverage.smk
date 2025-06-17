# Author: Dhatri Badri

rule coverage:
    input:
        r1 = lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R1.fastq.gz"),
        r2 = lambda wc: os.path.join(SHORT_READS_DIR, f"{wc.sample}_R2.fastq.gz"),
    output:
        coverage = "results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json",
    params:
        size = config["genome_size"]
    singularity:
        "docker://staphb/fastq-scan:1.0.1"
    shell:
        "zcat {input.r1} {input.r2} | fastq-scan -g {params.size} > {output.coverage}"