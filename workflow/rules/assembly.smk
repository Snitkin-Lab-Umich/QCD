# Author: Dhatri Badri

rule assembly:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        spades_assembly = f"results/{{prefix}}/spades/{{sample}}/contigs.fasta", 
        #updated_samples_file_successsfully = f"results/{{prefix}}/sample_files/updated_samples_passed_assembly.done",
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        #db = config["kraken_db"],
        threads= 10,
        prefix_dir = "results/{prefix}/",
        sample = "{sample}"
    #conda:
    #    "envs/spades.yaml"
    # singularity:
    #     "docker://staphb/spades:4.0.0"
    #envmodules:
    #    "Bioinformatics",
    #    "spades/4.0.0"
    # shell:
    #     "spades.py --isolate --pe1-1 {input.r1} --pe1-2 {input.r2} -o {params.out_dir} --threads {params.threads}" 
    wrapper:
        "file:workflow/wrapper_functions/spades"