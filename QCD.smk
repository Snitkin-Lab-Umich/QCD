# Author: Ali Pirani and Dhatri Badri
configfile: "config/config.yaml"

include: "QCD_report.smk"

import pandas as pd
import os

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

PREFIX = config["prefix"]

SHORTREADS = list(samples_df['sample_id'])

if not os.path.exists("results/"):
    os.system("mkdir %s" % "results/")

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)

def downsample_reads(file, file2, out1, out2, genome_size):
    file = file.pop()
    file2 = file2.pop()
    out1 = out1.pop()
    out2 = out2.pop()

    # Extract basic fastq reads stats with seqtk

    gsize = genome_size.pop()

    print("Using Genome Size: %s to calculate coverage" % gsize)
    
    seqtk_check = "/nfs/esnitkin/bin_group/seqtk/seqtk fqchk -q3 %s > %s_fastqchk.txt" % (file, file)

    print(seqtk_check)

    try:
        os.system(seqtk_check)
    except sp.CalledProcessError:
        print('Error running seqtk for extracting fastq statistics.')
        sys.exit(1)

    with open("%s_fastqchk.txt" % file, 'rU') as file_open:
        for line in file_open:
            if line.startswith('min_len'):
                line_split = line.split(';')
                min_len = line_split[0].split(': ')[1]
                max_len = line_split[1].split(': ')[1]
                avg_len = line_split[2].split(': ')[1]
            if line.startswith('ALL'):
                line_split = line.split('\t')
                total_bases = int(line_split[1]) * 2
    file_open.close()

    print('Average Read Length: %s' % avg_len)

    print('Total number of bases in fastq: %s' % total_bases)

    # Calculate original depth and check if it needs to be downsampled to a default coverage.
    ori_coverage_depth = int(total_bases / gsize)

    print('Original Covarage Depth: %s x' % ori_coverage_depth)

    # proc = sp.Popen(["nproc"], stdout=sp.PIPE, shell=True)
    # (nproc, err) = proc.communicate()
    # nproc = nproc.strip()

    if ori_coverage_depth >= 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        # r1_sub = "/tmp/%s" % os.path.basename(file)
        r1_sub = out1

        # Downsample using seqtk
        try:
            #print("Generating seqtk Downsampling command")
            print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file, factor, r1_sub))

            seqtk_downsample = "/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (
                file, factor, r1_sub)
            os.system(seqtk_downsample)
            #call(seqtk_downsample, logger)
        except sp.CalledProcessError:
            print('Error running seqtk for downsampling raw fastq reads.')
            sys.exit(1)

        if file2:
            # r2_sub = "/tmp/%s" % os.path.basename(file2)
            r2_sub = out2
            try:
                print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file2, factor, r2_sub))
                os.system("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (file2, factor, r2_sub))
                #call("seqtk sample %s %s | pigz --fast -c -p %s > /tmp/%s" % (file2, factor, nproc, os.path.basename(file2)), logger)
            except sp.CalledProcessError:
                print('Error running seqtk for downsampling raw fastq reads.')
                sys.exit(1)
        else:
            r2_sub = "None"

    elif ori_coverage_depth < 100:
        r1_sub = file
        r2_sub = file2
        os.system("cp %s %s" % (file, out1))
        os.system("cp %s %s" % (file2, out2))

rule all:
    input:
        coverage = expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", sample=SAMPLE, prefix=PREFIX),
        fastqc_raw = expand("results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html", sample=SAMPLE, prefix=PREFIX),
        trim = expand("results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE, prefix=PREFIX),
        fastqc_aftertrim = expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html", sample=SAMPLE, prefix=PREFIX),
        downsample_read = expand("results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE, prefix=PREFIX),
        spades_assembly = expand("results/{prefix}/spades/{sample}/contigs.fasta", sample=SAMPLE, prefix=PREFIX),
        spades_l1000_assembly = expand("results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta", sample=SAMPLE, prefix=PREFIX),
        prokka_gff = expand("results/{prefix}/prokka/{sample}/{sample}.gff", sample=SAMPLE, prefix=PREFIX),
        quast_report = expand("results/{prefix}/quast/{sample}/report.tsv", sample=SAMPLE, prefix=PREFIX),
        mlst_report = expand("results/{prefix}/mlst/{sample}/report.tsv", sample=SAMPLE, prefix=PREFIX),
        skani_ref_genome_results = expand("results/{prefix}/skani/{sample}/{sample}_skani_output.txt", sample=SAMPLE, prefix=PREFIX),
        coverage_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Final_Coverage.txt", prefix=PREFIX),
        skani_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Skani_report_final.csv", prefix=PREFIX),
        multiqc_report = expand("results/{prefix}/{prefix}_Report/multiqc/{prefix}_QC_report.html", prefix=PREFIX),
        mlst_final_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_MLST_results.csv", prefix=PREFIX),
        QC_summary = expand("results/{prefix}/{prefix}_Report/data/{prefix}_QC_summary.csv", prefix=PREFIX),
        QC_plot = expand("results/{prefix}/{prefix}_Report/fig/{prefix}_Coverage_distribution.png", prefix=PREFIX)
        
rule coverage:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
    output:
        coverage = f"results/{{prefix}}/raw_coverage/{{sample}}/{{sample}}_coverage.json",
    params:
        size=config["genome_size"]
    #conda:
    #    "envs/fastq-scan.yaml"
    singularity:
        "docker://staphb/fastq-scan:1.0.1"
    shell:
        "zcat {input.r1} {input.r2} | fastq-scan -g {params.size} > {output.coverage}"

rule quality_raw:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
    output:
        raw_fastqc_report_fwd = f"results/{{prefix}}/quality_raw/{{sample}}/{{sample}}_Forward/{{sample}}_R1_fastqc.html",
        raw_fastqc_report_rev = f"results/{{prefix}}/quality_raw/{{sample}}/{{sample}}_Reverse/{{sample}}_R2_fastqc.html",
    log:
        "logs/{prefix}/quality_raw/{sample}/{sample}.log"
    params:
        outdir="results/{prefix}/quality_raw/{sample}/{sample}"
    #conda:
    #    "envs/fastqc.yaml"
    singularity:
        "docker://staphb/fastqc:0.12.1"
    #envmodules:
    #    "Bioinformatics",
    #    "fastqc"
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        """

rule trimmomatic_pe:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz"))  
    output:
        r1 = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R1_trim_paired.fastq.gz",
        r2 = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R2_trim_paired.fastq.gz", 
        # reads where trimming entirely removed the mate
        r1_unpaired = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = f"results/{{prefix}}/trimmomatic/{{sample}}/{{sample}}_R2_trim_unpaired.fastq.gz",
    params:
        adapter_filepath=config["adapter_file"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],
        threads = config["ncores"],
    log:
        "logs/{prefix}/trimmomatic/{sample}/{sample}.log"
    #conda:
    #    "envs/trimmomatic.yaml"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    #envmodules:
    #    "Bioinformatics",
    #    "trimmomatic"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {params.threads} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}"

rule quality_aftertrim:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        aftertrim_fastqc_report_fwd = f"results/{{prefix}}/quality_aftertrim/{{sample}}/{{sample}}_Forward/{{sample}}_R1_trim_paired_fastqc.html",
        aftertrim_fastqc_report_rev = f"results/{{prefix}}/quality_aftertrim/{{sample}}/{{sample}}_Reverse/{{sample}}_R2_trim_paired_fastqc.html",
    log:
        "logs/{prefix}/{sample}/quality_aftertrim/{sample}.log"
    params:
        outdir="results/{prefix}/quality_aftertrim/{sample}/{sample}"
    #conda:
    #    "envs/fastqc.yaml"
    singularity:
        "docker://staphb/fastqc:0.12.1"
    #envmodules:
    #    "Bioinformatics",
    #    "fastqc"
    shell:
        """
        mkdir -p {params.outdir}_Forward {params.outdir}_Reverse
        fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}
        """

rule downsample:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        outr1 = f"results/{{prefix}}/downsample/{{sample}}/{{sample}}_R1_trim_paired.fastq.gz",
        outr2 = f"results/{{prefix}}/downsample/{{sample}}/{{sample}}_R2_trim_paired.fastq.gz",
    params:
        gsize = config["genome_size"],
    run:
        downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize})

rule kraken:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        kraken_out = f"results/{{prefix}}/kraken/{{sample}}/{{sample}}_kraken_out",
        kraken_report = f"results/{{prefix}}/kraken/{{sample}}/{{sample}}_kraken_report.tsv",
    params:
        db = config["kraken_db"],
        threads = 12
        # threads = config["threads"]
    #conda:
    #    "envs/kraken.yaml"
    singularity:
        "docker://staphb/kraken2:2.1.3"
    shell:
        "kraken2 --db {params.db} --threads {params.threads} --paired --gzip-compressed --quick --output {output.kraken_out} --report {output.kraken_report} {input.r1} {input.r2}"

rule assembly:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/downsample/{wildcards.sample}/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        spades_assembly = f"results/{{prefix}}/spades/{{sample}}/contigs.fasta",
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        db = config["kraken_db"],
    #conda:
    #    "envs/spades.yaml"
    singularity:
        "docker://staphb/spades:4.0.0"
    #envmodules:
    #    "Bioinformatics",
    #    "spades/4.0.0"
    shell:
        "spades.py --isolate --pe1-1 {input.r1} --pe1-2 {input.r2} -o {params.out_dir}"

#rule bioawk:
#    input:
#        spades_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/contigs.fasta"),
#    output:
#        spades_l1000_assembly = f"results/{{prefix}}/{{sample}}/spades/{{sample}}_contigs_l1000.fasta",
#    params:
#        out_dir = "results/{prefix}/{sample}/spades/",
#        bioawk_params = config["bioawk"],
#        prefix = "{sample}",
#    conda:
#        "envs/bioawk.yaml"
#    shell:
#        """
#        {params.bioawk_params} {input.spades_assembly} > {output.spades_l1000_assembly} && grep '>' {output.spades_l1000_assembly} > {params.out_dir}/spades_assembly_header_info.txt && sed -i 's/>NODE_/>{params.prefix}_/g' {output.spades_l1000_assembly} && sed -i 's/_length_.*_cov_.*//g' {output.spades_l1000_assembly}
#        """

rule bioawk:
    input:
        spades_assembly = lambda wildcards: f"results/{wildcards.prefix}/spades/{wildcards.sample}/contigs.fasta"
    output:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta"
    params:
        out_dir = "results/{prefix}/spades/{sample}/",
        bioawk_params = config["bioawk"],
        prefix = "{sample}"
    conda:
        "envs/bioawk.yaml"
    shell:
        """
        ./bioawk.sh {input.spades_assembly} {output.spades_l1000_assembly} {params.out_dir} {params.prefix}
        """

rule prokka:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        prokka_gff = f"results/{{prefix}}/prokka/{{sample}}/{{sample}}.gff",
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

rule quast:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        quast_report = f"results/{{prefix}}/quast/{{sample}}/report.tsv",
    params: 
        outdir = "results/{prefix}/quast/{sample}/",
        prefix = "{sample}",
    #conda:
    #    "envs/quast.yaml"
    singularity:
        "docker://staphb/quast:5.0.2"
    #envmodules:
    #    "Bioinformatics",
    #    "quast"
    shell:
        """
        ./quast.sh {input.spades_l1000_assembly} {params.outdir} 
        """
        #"quast.py {input.spades_l1000_assembly} -o {params.outdir} --contig-thresholds 0,1000,5000,10000,25000,50000"

rule mlst:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        mlst_report = f"results/{{prefix}}/mlst/{{sample}}/report.tsv",
    params: 
        outdir = "results/{prefix}/mlst/{sample}/",
        prefix = "{sample}",
    #conda:
    #    "envs/mlst.yaml"
    singularity:
        "docker://staphb/mlst:2.23.0-2024-03"
    #envmodules:
    #    "Bioinformatics",
    #    "mlst"
    shell:
        "mlst {input.spades_l1000_assembly} > {output.mlst_report}"

#rule amrfinder:
#    input:
#        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/{wildcards.sample}_contigs_l1000.fasta"),
#    output:
#        amrfinder = f"results/{{prefix}}/{{sample}}/amrfinder/{{sample}}_amrfinder.tsv",
#    params: 
#        outdir = "results/{prefix}/{sample}/amrfinder",
#        prefix = "{sample}",
#        organism = config['amrfinder_organism']
    #conda:
    #    "envs/amrfinder.yaml"
#    singularity:
#        "docker://staphb/ncbi-amrfinderplus:latest"
#    shell:
#        "amrfinder --plus --output {output.amrfinder} --debug --log {params.outdir}/{params.prefix}.log --nucleotide_output {params.outdir}/{params.prefix}_reported_nucl.fna -n {input.spades_l1000_assembly} -O {params.organism}"

rule busco:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        busco_out = f"results/{{prefix}}/{{sample}}/busco/busco.txt",
    params: 
        outdir = "results/{prefix}/busco/{sample}/",
        prefix = "{sample}",
        threads = config["ncores"],
    #conda:
    #    "envs/busco.yaml"
    singularity:
        "docker://staphb/busco:5.7.1-prok-bacteria_odb10_2024-01-08"
    #envmodules:
    #    "Bioinformatics",
    #    "busco"
    shell:
        "busco -f -i {input.spades_l1000_assembly} -m genome -l bacteria_odb10 -o {params.outdir}"

rule skani:
    input:
        spades_contigs_file = lambda wildcards: expand(f"results/{wildcards.prefix}/spades/{wildcards.sample}/contigs.fasta")
    output:
        skani_output = f"results/{{prefix}}/skani/{{sample}}/{{sample}}_skani_output.txt"
    params:
        skani_ani_db = config["skani_db"],
        threads = 4
    #conda:
    #    "envs/skani.yaml"
    singularity:
        "docker://staphb/skani:0.2.1"
    shell:
        "skani search {input.spades_contigs_file} -d {params.skani_ani_db} -o {output.skani_output} -t {params.threads}"
        
#rule multiqc:
#    input:
#        quast_report = f"results/{{prefix}}/{{sample}}/quast/report.tsv",
#        prokka_gff = f"results/{{prefix}}/{{sample}}/prokka/{{sample}}.gff",
#        spades_assembly = f"results/{{prefix}}/{{sample}}/spades/contigs.fasta",
#        kraken_report = f"results/{{prefix}}/{{sample}}/kraken/{{sample}}_kraken_report.tsv",
#        coverage = f"results/{{prefix}}/{{sample}}/raw_coverage/{{sample}}_coverage.json",
#        aftertrim_fastqc_report_fwd = f"results/{{prefix}}/{{sample}}/quality_aftertrim/{{sample}}_Forward/{{sample}}_R1_trim_paired_fastqc.html",
#        aftertrim_fastqc_report_rev = f"results/{{prefix}}/{{sample}}/quality_aftertrim/{{sample}}_Reverse/{{sample}}_R2_trim_paired_fastqc.html",
#        raw_fastqc_report_fwd = f"results/{{prefix}}/{{sample}}/quality_raw/{{sample}}_Forward/{{sample}}_R1_fastqc.html",
#        raw_fastqc_report_rev = f"results/{{prefix}}/{{sample}}/quality_raw/{{sample}}_Reverse/{{sample}}_R2_fastqc.html"
#    output:
#        multiqc_report = f"results/{{prefix}}/multiqc/{{prefix}}_QC_report.html",
#    params:
#        resultsoutdir = "results/{prefix}",
#        outdir = "results/{prefix}/multiqc",
#        prefix = "{prefix}",
    #conda:
    #    "envs/multiqc.yaml"
    #singularity:
        #"docker://staphb/multiqc:1.18"
#    shell:
#        """
#        module load Bioinformatics
#        module load multiqc
#        multiqc -f --outdir {params.outdir} -n {params.prefix}_QC_report -i {params.prefix}_QC_report {params.resultsoutdir}
#        """
        
