rule busco:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        busco_out = "results/{prefix}/busco/{sample}/{sample}_busco_out.txt",
    params: 
        outdir = "results/{prefix}/busco/{sample}/",
        prefix = "{sample}",
        threads = config["ncores"],
        assembly_busco_out = "short_summary.specific.bacteria_odb10.{sample}.txt"
    #conda:
    #    "envs/busco.yaml"
    singularity:
        "docker://staphb/busco:6.0.0-prok-bacteria_odb12_2024-11-14"
    #envmodules:
    #    "Bioinformatics",
    #    "busco"
    shell:
        """
        busco -f -i {input.spades_l1000_assembly} -m genome -l bacteria_odb10 -o {params.outdir}
        cp {params.outdir}/{params.assembly_busco_out} {params.outdir}/{params.prefix}_busco_out.txt
        """
