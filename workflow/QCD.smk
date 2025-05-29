# Author: Ali Pirani and Dhatri Badri
configfile: "config/config.yaml"

import pandas as pd
import os
import json

# Run workflow until coverage
samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])
PREFIX = config["prefix"]

if not os.path.exists(f"results/{PREFIX}/sample_files"):
    os.makedirs(f"results/{PREFIX}/sample_files")

# Run workflow until assembly
samps_passed_cov = f"results/{PREFIX}/sample_files/samples_passed_coverage.csv"
if not os.path.exists(samps_passed_cov):
    with open(samps_passed_cov, 'w') as spc:
        _ = spc.write('sample_id,illumina_r1\n')

samples_df_passed_cov = pd.read_csv(samps_passed_cov)
SAMPLE_PASS_COV = list(samples_df_passed_cov['sample_id'])

# Run workflow post assembly
samps_passed_assembly = f"results/{PREFIX}/sample_files/samples_passed_assembly.csv"
if not os.path.exists(samps_passed_assembly):
    with open(samps_passed_assembly, 'w') as spa:
        _ = spa.write('sample_id,illumina_r1\n')

samples_df_passed_assembly = pd.read_csv(samps_passed_assembly) 
SAMPLE_PASS_ASSEMBLY = list(samples_df_passed_assembly['sample_id'])



rule all:
    input:
        coverage= expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", sample=SAMPLE, prefix=PREFIX),
        samples_passed_coverage = expand("results/{prefix}/sample_files/samples_passed_coverage.csv", prefix=PREFIX), 
        updated_samples_passed_coverage = expand("results/{prefix}/sample_files/updated_samples_passed_coverage.done", prefix=PREFIX),
        ####
        fastqc_raw = expand("results/{prefix}/quality_raw/{sample}/{sample}_Forward/{sample}_R1_fastqc.html", sample=SAMPLE_PASS_COV, prefix=PREFIX),
        trim = expand("results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE_PASS_COV, prefix=PREFIX),
        fastqc_aftertrim = expand("results/{prefix}/quality_aftertrim/{sample}/{sample}_Forward/{sample}_R1_trim_paired_fastqc.html", sample=SAMPLE_PASS_COV, prefix=PREFIX),
        downsample_read = expand("results/{prefix}/downsample/{sample}/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE_PASS_COV, prefix=PREFIX),
        spades_assembly = expand("results/{prefix}/spades/{sample}/contigs.fasta", sample=SAMPLE_PASS_COV, prefix=PREFIX),
        samples_passed_assembly = expand("results/{prefix}/sample_files/samples_passed_assembly.csv", prefix=PREFIX),
        #####
        spades_l1000_assembly = expand("results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta", sample=SAMPLE_PASS_ASSEMBLY, prefix=PREFIX),
        prokka_gff = expand("results/{prefix}/prokka/{sample}/{sample}.gff", sample=SAMPLE_PASS_ASSEMBLY, prefix=PREFIX),
        quast_report = expand("results/{prefix}/quast/{sample}/report.tsv", sample=SAMPLE_PASS_ASSEMBLY, prefix=PREFIX),
        sample_mlst_report = expand("results/{prefix}/mlst/{sample}/report.tsv", sample=SAMPLE_PASS_ASSEMBLY, prefix=PREFIX),
        skani_ref_genome_results = expand("results/{prefix}/skani/{sample}/{sample}_skani_output.txt", sample=SAMPLE_PASS_ASSEMBLY, prefix=PREFIX),

################################################################################################################################################################################################
rule coverage:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
    output:
        coverage = "results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json",
    params:
        size = config["genome_size"]
    singularity:
        "docker://staphb/fastq-scan:1.0.1"
    shell:
        "zcat {input.r1} {input.r2} | fastq-scan -g {params.size} > {output.coverage}"

def samples_that_passed_coverage(coverage_json, passed_csv, updated_samples_file, outdir):
    with open(coverage_json) as f:
        data = json.load(f)

    sample_name = os.path.basename(coverage_json).replace("_coverage.json", "")
    failed_csv = os.path.join(outdir, "sample_files/samples_failed_coverage_summary.csv")

    coverage = data.get('qc_stats', {}).get('coverage', 0)
    total_reads = data.get('qc_stats', {}).get('read_total', 0)
    total_bp = data.get('qc_stats', {}).get('total_bp', 0)
    mean_read_length = data.get('qc_stats', {}).get('read_mean', 0)

    # Ensure passed_csv exists before appending
    if not os.path.exists((passed_csv)):
        pd.DataFrame(columns=["sample_id", "illumina_r1"]).to_csv(passed_csv, index=False)

    # Read existing passed samples to avoid duplicates
    existing_passed_samples = set(pd.read_csv(passed_csv)["sample_id"])

    if coverage > 20 and sample_name not in existing_passed_samples:
        passed_df = pd.DataFrame([[sample_name, f"{sample_name}_R1.fastq.gz"]],
                                 columns=["sample_id", "illumina_r1"])
        passed_df.to_csv((passed_csv), mode='a', index=False, header=False)

        with open(updated_samples_file, "a") as done_file:
            done_file.write(f"Sample: {sample_name} updated successfully.\n")

    elif coverage <= 20:
        failed_df = pd.DataFrame([[sample_name, total_reads, total_bp, mean_read_length, coverage] + ["NA"] * 15 + ["FAIL"] + ["NA"] * 5],
                                 columns=[
                                     "Sample", "Total_reads", "Total_bp", "MeanReadLength", "Coverage", "Scheme",
                                     "ST", "After_trim_per_base_sequence_content", "After_trim_overrepresented_sequences",
                                     "After_trim_%GC", "After_trim_Total Bases", "After_trim_Total Sequences",
                                     "After_trim_median_sequence_length", "After_trim_avg_sequence_length",
                                     "After_trim_total_deduplicated_percentage", "After_trim_Sequence length",
                                     "After_trim_adapter_content", "N50", "Total length", "Total # of contigs",
                                     "QC Check", "ANI", "Align_fraction_ref", "Align_fraction_query",
                                     "Ref_name", "Species"
                                 ])
        failed_df.to_csv(failed_csv, mode='a', index=False, header=not os.path.exists(failed_csv))

rule create_samples_passed_coverage_file:
    input:
        coverage = expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", prefix=PREFIX, sample=SAMPLE),
        samples_passed_coverage = expand("results/{prefix}/sample_files/samples_passed_coverage.csv", prefix=PREFIX,),
    output:
        #samples_passed_coverage = "results/{prefix}/sample_files/samples_passed_coverage.csv",
        updated_samples_file_successsfully = "results/{prefix}/sample_files/updated_samples_passed_coverage.done"
    params:
        outdir = "results/{prefix}"
    run:
        for coverage_json in input.coverage:
            samples_that_passed_coverage(coverage_json, input.samples_passed_coverage[0], output.updated_samples_file_successsfully, params.outdir)

#####################################################################################################################################

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
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
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
    # run:
    #     downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize})
    wrapper:
        "file:workflow/wrapper_functions/downsample"

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
 
#####################################################################################################################################################################

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

rule quast:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        quast_report = "results/{prefix}/quast/{sample}/report.tsv",
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
        workflow/scripts/quast.sh {input.spades_l1000_assembly} {params.outdir} 
        """
        #"quast.py {input.spades_l1000_assembly} -o {params.outdir} --contig-thresholds 0,1000,5000,10000,25000,50000"

rule mlst:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        mlst_report = "results/{prefix}/mlst/{sample}/report.tsv",
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


rule busco:
    input:
        spades_l1000_assembly = "results/{prefix}/spades/{sample}/{sample}_contigs_l1000.fasta",
    output:
        busco_out = "results/{prefix}/{sample}/busco/busco.txt",
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

#######################################################################################################################################################
"""
END OF PIPELINE
"""