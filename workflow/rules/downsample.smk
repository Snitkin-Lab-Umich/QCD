# Author: Dhatri Badri

rule downsample:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        outr1 = "results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz",
        outr2 = "results/{prefix}/downsample/{sample}/{sample}_R2_trim_paired.fastq.gz",
    params:
        gsize = config["genome_size"],
    # run:
    #     downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize})
    wrapper:
        "file:workflow/wrapper_functions/downsample"