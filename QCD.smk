# Author: Ali Pirani
configfile: "config/config.yaml"

import pandas as pd
import os

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

PREFIX = config["prefix"]

SHORTREADS = list(samples_df['sample_id'])

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
        

    # # return r1_sub, r2_sub, seqtk_downsample


rule all:
    input:
        coverage = expand("results/{prefix}/{sample}/raw_coverage/{sample}_coverage.json", sample=SAMPLE, prefix=PREFIX),
        fastqc_raw = expand("results/{prefix}/{sample}/quality_raw/{sample}_Forward/{sample}_R1_fastqc.html", sample=SAMPLE, prefix=PREFIX),
        trim = expand("results/{prefix}/{sample}/trimmomatic/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE, prefix=PREFIX),
        fastqc_aftertrim = expand("results/{prefix}/{sample}/quality_aftertrim/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html", sample=SAMPLE, prefix=PREFIX),
        #ariba_report = expand("results/{prefix}/{sample}/ariba_card/report.tsv", sample=SAMPLE, prefix=PREFIX),
        #ariba_mlst_report = expand("results/{prefix}/{sample}/ariba_mlst/mlst_report.tsv", sample=SAMPLE, prefix=PREFIX),
        downsample_read = expand("results/{prefix}/{sample}/downsample/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE, prefix=PREFIX),
        spades_assembly = expand("results/{prefix}/{sample}/spades/contigs.fasta", sample=SAMPLE, prefix=PREFIX),
        spades_l1000_assembly = expand("results/{prefix}/{sample}/spades/{sample}_contigs_l1000.fasta", sample=SAMPLE, prefix=PREFIX),
        prokka_gff = expand("results/{prefix}/{sample}/prokka/{sample}.gff", sample=SAMPLE, prefix=PREFIX),
        quast_report = expand("results/{prefix}/{sample}/quast/report.tsv", sample=SAMPLE, prefix=PREFIX),
        amrfinder = expand("results/{prefix}/{sample}/amrfinder/{sample}_amrfinder.tsv", sample=SAMPLE, prefix=PREFIX),
        mlst_report = expand("results/{prefix}/{sample}/mlst/report.tsv", sample=SAMPLE, prefix=PREFIX),
        kraken_report = expand("results/{prefix}/{sample}/kraken/{sample}_kraken_report.tsv", sample=SAMPLE, prefix=PREFIX),

rule coverage:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
    output:
        coverage = f"results/{{prefix}}/{{sample}}/raw_coverage/{{sample}}_coverage.json",
    params:
        size=config["genome_size"]
    conda:
        "envs/fastq-scan.yaml"
    shell:
        "zcat {input.r1} {input.r2} | fastq-scan -g {params.size} > {output.coverage}"

rule quality_raw:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
    output:
        raw_fastqc_report_fwd = f"results/{{prefix}}/{{sample}}/quality_raw/{{sample}}_Forward/{{sample}}_R1_fastqc.html",
        raw_fastqc_report_rev = f"results/{{prefix}}/{{sample}}/quality_raw/{{sample}}_Reverse/{{sample}}_R2_fastqc.html",
    log:
        "logs/{prefix}/{sample}/quality_raw/{sample}.log"
    params:
        outdir="results/{prefix}/{sample}/quality_raw/{sample}"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}"

rule trimmomatic_pe:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        
    output:
        r1 = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_paired.fastq.gz",
        r2 = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_paired.fastq.gz", 
        # reads where trimming entirely removed the mate
        r1_unpaired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_unpaired.fastq.gz",
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
        "logs/{prefix}/{sample}/trimmomatic/{sample}.log"
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {params.threads} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}"

rule quality_aftertrim:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        aftertrim_fastqc_report_fwd = f"results/{{prefix}}/{{sample}}/quality_aftertrim/{{sample}}_Forward/{{sample}}_R1_trim_paired_fastqc.html",
        aftertrim_fastqc_report_rev = f"results/{{prefix}}/{{sample}}/quality_aftertrim/{{sample}}_Reverse/{{sample}}_R2_trim_paired_fastqc.html",
    log:
        "logs/{prefix}/{sample}/quality_aftertrim/{sample}.log"
    params:
        outdir="results/{prefix}/{sample}/quality_aftertrim/{sample}"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc -o {params.outdir}_Forward {input.r1} && fastqc -o {params.outdir}_Reverse {input.r2} &>{log}"

rule ariba_card:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),

    output:
        ariba_report = f"results/{{prefix}}/{{sample}}/ariba_card/report.tsv",
    params:
        card_db = config["ariba_card_db"],
        threads = config["ncores"],
        outdir = "results/{prefix}/{sample}/ariba_card/"
    conda:
        "envs/ariba.yaml"
    shell:
        "ariba run --verbose --force {params.card_db} {input.r1} {input.r2} {params.outdir} --tmp_dir /tmp/"

rule ariba_mlst:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),

    output:
        ariba_mlstreport = f"results/{{prefix}}/{{sample}}/ariba_mlst/mlst_report.tsv",
    params:
        mlst_db = config["ariba_mlst_db"],
        threads = config["ncores"],
        outdir = "results/{prefix}/{sample}/ariba_mlst/"
    conda:
        "envs/ariba.yaml"
    shell:
        "ariba run --verbose --force {params.mlst_db} {input.r1} {input.r2} {params.outdir} --tmp_dir /tmp/"

rule ariba_summary:
    input:
        ariba_cardreport = lambda wildcards: expand(f"results/{wildcards.prefix}/*/ariba_card/report.tsv"),
    output:
        ariba_cardsummary = f"results/{{prefix}}/Ariba_Card_Minimal_Summary.csv",
    
    params:
        prefix = "results/{prefix}/{prefix}_Ariba_Card_Minimal_Summary",
        threads = config["ncores"],
    conda:
        "envs/ariba.yaml"
    shell:
        "ariba summary {params.prefix} {input.ariba_cardreport} > {output.ariba_cardsummary}"

rule downsample:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),

    output:
        outr1 = f"results/{{prefix}}/{{sample}}/downsample/{{sample}}_R1_trim_paired.fastq.gz",
        outr2 = f"results/{{prefix}}/{{sample}}/downsample/{{sample}}_R2_trim_paired.fastq.gz",
    params:
        gsize = config["genome_size"],
    run:
        downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize})

rule kraken:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/downsample/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/downsample/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        kraken_out = f"results/{{prefix}}/{{sample}}/kraken/{{sample}}_kraken_out",
        kraken_report = f"results/{{prefix}}/{{sample}}/kraken/{{sample}}_kraken_report.tsv",
    params:
        db = config["kraken_db"],
        threads = 12
        # threads = config["threads"]
    conda:
        "envs/kraken.yaml"
    shell:
        "kraken2 --db {params.db} --threads {params.threads} --paired --gzip-compressed --quick --output {output.kraken_out} --report {output.kraken_report} {input.r1} {input.r2}"

rule assembly:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/downsample/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/downsample/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        spades_assembly = f"results/{{prefix}}/{{sample}}/spades/contigs.fasta",
    params:
        out_dir = "results/{prefix}/{sample}/spades/",
        db = config["kraken_db"],
    conda:
        "envs/spades.yaml"
    shell:
        "spades.py --careful --sc --pe1-1 {input.r1} --pe1-2 {input.r2} -o {params.out_dir}"

rule bioawk:
    input:
        spades_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/contigs.fasta"),
    output:
        spades_l1000_assembly = f"results/{{prefix}}/{{sample}}/spades/{{sample}}_contigs_l1000.fasta",
    params:
        out_dir = "results/{prefix}/{sample}/spades/",
        bioawk_params = config["bioawk"],
        prefix = "{sample}",
    conda:
        "envs/bioawk.yaml"
    shell:
        "{params.bioawk_params} {input.spades_assembly} > {output.spades_l1000_assembly} && grep '>' {output.spades_l1000_assembly} > {params.out_dir}/spades_assembly_header_info.txt && sed -i 's/>NODE_/>{params.prefix}_/g' {output.spades_l1000_assembly} && sed -i 's/_length_.*_cov_.*//g' {output.spades_l1000_assembly}"

rule prokka:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        prokka_gff = f"results/{{prefix}}/{{sample}}/prokka/{{sample}}.gff",
    params: 
        prokka_params = config["prokka"],
        outdir = "results/{prefix}/{sample}/prokka",
        prefix = "{sample}",
    conda:
        "envs/prokka.yaml"
    shell:
        "prokka -outdir {params.outdir} --strain {params.prefix} --prefix {params.prefix} {params.prokka_params} {input.spades_l1000_assembly}"


rule quast:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        quast_report = f"results/{{prefix}}/{{sample}}/quast/report.tsv",
    params: 
        outdir = "results/{prefix}/{sample}/quast",
        prefix = "{sample}",
    conda:
        "envs/quast.yaml"
    shell:
        "quast.py {input.spades_l1000_assembly} -o {params.outdir} --contig-thresholds 0,1000,5000,10000,25000,50000"

rule mlst:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        mlst_report = f"results/{{prefix}}/{{sample}}/mlst/report.tsv",
    params: 
        outdir = "results/{prefix}/{sample}/mlst",
        prefix = "{sample}",
    conda:
        "envs/mlst.yaml"
    shell:
        "mlst {input.spades_l1000_assembly} > {output.mlst_report}"

rule amrfinder:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        amrfinder = f"results/{{prefix}}/{{sample}}/amrfinder/{{sample}}_amrfinder.tsv",
    params: 
        outdir = "results/{prefix}/{sample}/amrfinder",
        prefix = "{sample}",
        organism = config['amrfinder_organism']
    conda:
        "envs/amrfinder.yaml"
    shell:
        "amrfinder --plus --output {output.amrfinder} --debug --log {params.outdir}/{params.prefix}.log --nucleotide_output {params.outdir}/{params.prefix}_reported_nucl.fna -n {input.spades_l1000_assembly}"
        

rule busco:
    input:
        spades_l1000_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/spades/{wildcards.sample}_contigs_l1000.fasta"),
    output:
        busco_out = f"results/{{prefix}}/{{sample}}/busco/busco.txt",
    params: 
        outdir = "results/{prefix}/{sample}/busco",
        prefix = "{sample}",
        threads = config["ncores"],
    conda:
        "envs/busco.yaml"
    shell:
        "busco -f -i {input.spades_l1000_assembly} -m genome -l bacteria_odb10 -o {params.outdir}"
        
rule multiqc:
    input:
        quast_report = f"results/{{prefix}}/{{sample}}/quast/report.tsv",
        prokka_gff = f"results/{{prefix}}/{{sample}}/prokka/{{sample}}.gff",
        spades_assembly = f"results/{{prefix}}/{{sample}}/spades/contigs.fasta",
        kraken_report = f"results/{{prefix}}/{{sample}}/kraken/{{sample}}_kraken_report.tsv",
        coverage = f"results/{{prefix}}/{{sample}}/raw_coverage/{{sample}}_coverage.json",
        aftertrim_fastqc_report_fwd = f"results/{{prefix}}/{{sample}}/quality_aftertrim/{{sample}}_Forward/{{sample}}_R1_trim_paired_fastqc.html",
        aftertrim_fastqc_report_rev = f"results/{{prefix}}/{{sample}}/quality_aftertrim/{{sample}}_Reverse/{{sample}}_R2_trim_paired_fastqc.html",
        raw_fastqc_report_fwd = f"results/{{prefix}}/{{sample}}/quality_raw/{{sample}}_Forward/{{sample}}_R1_fastqc.html",
        raw_fastqc_report_rev = f"results/{{prefix}}/{{sample}}/quality_raw/{{sample}}_Reverse/{{sample}}_R2_fastqc.html",
        
    output:
        multiqc_report = f"results/{{prefix}}/multiqc/{{prefix}}_QC_report.html",
    params:
        resultsoutdir = "results/{prefix}",
        outdir = "results/{prefix}/multiqc",
        prefix = "{prefix}",
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc -f --outdir {params.outdir} -n {params.prefix}_QC_report -i {params.prefix}_QC_report {params.resultsoutdir}"
        
