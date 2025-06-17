rule bioawk:
    input:
        spades_assembly = "results/{prefix}/spades/{sample}/contigs.fasta"
    output:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta"
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        bioawk_params = config["bioawk"],
        prefix = "{sample}"
    conda:
        "envs/bioawk.yaml"
    # singularity:
    #     "docker://lbmc/bioawk:1.0"
    shell:
        """
        workflow/scripts/bioawk.sh {input.spades_assembly} {output.spades_l1000_assembly} {params.out_dir} {params.prefix}
        """