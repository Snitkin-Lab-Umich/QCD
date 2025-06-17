rule prokka:
    input:
        spades_l1000_assembly ="results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        prokka_gff = "results/{prefix}/prokka/{sample}/{sample}.gff",
    params: 
        prokka_params = config["prokka"],
        outdir = "results/{prefix}/prokka/{sample}",
        prefix = "{sample}",
    #conda:
    #    "envs/prokka.yaml"
    singularity:
        "docker://staphb/prokka:1.14.6"
    #envmodules:
    #    "Bioinformatics",
    #    "prokka"
    shell:
        "prokka -outdir {params.outdir} --strain {params.prefix} --prefix {params.prefix} {params.prokka_params} {input.spades_l1000_assembly}"
