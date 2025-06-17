rule skani:
    input:
        spades_contigs_file = "results/{prefix}/spades/{sample}/contigs.fasta"
    output:
        skani_output = "results/{prefix}/skani/{sample}/{sample}_skani_output.txt"
    params:
        skani_ani_db = config["skani_db"],
        threads = 6
    #conda:
    #    "envs/skani.yaml"
    singularity:
        "docker://staphb/skani:0.2.1"
    shell:
        "skani search {input.spades_contigs_file} -d {params.skani_ani_db} -o {output.skani_output} -t {params.threads}"
