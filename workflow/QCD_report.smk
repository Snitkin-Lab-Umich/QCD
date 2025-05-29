# Author: Ali Pirani and Dhatri Badri  
configfile: "config/config.yaml"

import pandas as pd
import os
import json
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re

PREFIX = config["prefix"]

samps_passed_assembly = f"results/{PREFIX}/sample_files/samples_passed_assembly.csv"

samples_df_passed_assembly = pd.read_csv(samps_passed_assembly)
SAMPLE_PASS_ASSEMBLY = list(samples_df_passed_assembly['sample_id'])

if not os.path.exists(f"results/{PREFIX}"):
    os.makedirs(f"results/{PREFIX}")

# Organize reports directory
prefix = PREFIX
outdir = "results/%s" % prefix
report_dir = outdir + "/%s_Report" % prefix
report_script_dir = report_dir + "/scripts"
report_data_dir = report_dir + "/data"
report_multiqc_dir = report_dir + "/multiqc"
report_fig_dir = report_dir + "/fig"

isExist = os.path.exists(report_dir)
if not isExist:
    os.makedirs(report_dir)

isExist = os.path.exists(report_script_dir)
if not isExist:
    os.makedirs(report_script_dir)

isExist = os.path.exists(report_data_dir)
if not isExist:
    os.makedirs(report_data_dir)

isExist = os.path.exists(report_multiqc_dir)
if not isExist:
    os.makedirs(report_multiqc_dir)

isExist = os.path.exists(report_fig_dir)
if not isExist:
    os.makedirs(report_fig_dir)

def coverage_report(prefix, outdir):
    prefix = prefix.pop()
    report_dir = str(outdir.pop()) + "/%s_Report" % prefix
    # Generate Coverage report 
    final_coverage_file = "%s/data/%s_Final_Coverage.txt" % (report_dir, prefix)
    f3=open(final_coverage_file, 'w+')
    header = "Sample,Total_reads,Total_bp,MeanReadLength,Coverage\n"
    f3.write(header)

    for sampl in SAMPLE_PASS_ASSEMBLY:
        coverage_json = "results/%s/raw_coverage/%s/%s_coverage.json" % (prefix, sampl, sampl)
        f = open(coverage_json)
        data = json.load(f)
        # data = json.loads(coverage_json)
        f3.write("%s,%s,%s,%s,%s\n" % (sampl, data['qc_stats']['read_total'], data['qc_stats']['total_bp'], data['qc_stats']['read_mean'], data['qc_stats']['coverage']))
    f3.close()   

    Coverage = pd.read_csv(final_coverage_file, sep=',', header=0)
    Coverage = Coverage.replace(['_R1.fastq.gz'], '', regex=True)

    #print ("Number of Samples in Coverage Report - %s" % len(Coverage))

    #Coverage_NEG_CNTL = Coverage[Coverage.Sample.str.match('(.*NEG*)')]

    #print ("Number of Negative Control samples %s" % len(Coverage_NEG_CNTL))

    #print ("Number of Negative Control samples with > 100X coverage %s" % len(Coverage_NEG_CNTL[Coverage_NEG_CNTL['Coverage'] > 100]))

    #Coverage_dist = Coverage.sort_values(by='Coverage',ascending=False).plot(x='Sample_name', y='Coverage', kind="barh", title="Estimated Genome Coverage", figsize=(20, 20), fontsize=40).get_figure()

    #Coverage_dist.savefig('%s/%s_Coverage_distribution.pdf' % (report_dir, prefix))


def skani_report(outdir, prefix):
    prefix = prefix.pop()
    outdir = "results/%s" % prefix
    report_dir = str(outdir) + "/%s_Report" % prefix
    report_data_dir = report_dir + "/data"
    result_df = pd.DataFrame(columns=['Sample', 'ANI', 'Align_fraction_ref', 'Align_fraction_query', 'Ref_name', 'Species'])  # Add 'Species' column

    skani_dir = os.path.join(outdir, 'skani')  # Navigate to skani directory

    for sample_name in os.listdir(skani_dir):  # Iterate over samples in the results/prefix/skani directory
        sample_dir = os.path.join(skani_dir, sample_name)

        if os.path.isdir(sample_dir):  # Check if it's a directory
            skani_file_path = os.path.join(sample_dir, f'{sample_name}_skani_output.txt')  # Look for the skani output file

            if os.path.exists(skani_file_path):  # Check if the skani file exists
                skani_file = pd.read_csv(skani_file_path, sep='\t| ,', skipinitialspace=True, header=0)  # Read the skani file
                first_row_df = skani_file[['ANI', 'Align_fraction_ref', 'Align_fraction_query', 'Ref_name']].iloc[:1]  # Extract the first row

                if first_row_df.empty:  # Check if the first row is empty
                    first_row_df = pd.DataFrame({
                        'Sample': [sample_name],  # Add sample name
                        'ANI': ["NA"],
                        'Align_fraction_ref': ["NA"],
                        'Align_fraction_query': ["NA"],
                        'Ref_name': ["NA"],
                        'Species': ["NA"]  # Add NAs for Species
                    })
                else:
                    first_row_df.loc[:, 'Sample'] = sample_name  # Add sample name
                    # Extract species using regex from Ref_name
                    first_row_df.loc[:, 'Species'] = first_row_df['Ref_name'].apply(
                        lambda x: re.search(r"[A-Za-z]+\s[A-Za-z]+", x).group(0) if pd.notnull(x) and re.search(r"[A-Za-z]+\s[A-Za-z]+", x) else "NAs"
                    )

                first_row_df = first_row_df[['Sample', 'ANI', 'Align_fraction_ref', 'Align_fraction_query', 'Ref_name', 'Species']]  # Reorder columns
                result_df = pd.concat([result_df, first_row_df], ignore_index=True)  # Concatenate to the result dataframe

    result_file_path = os.path.join(report_data_dir, f'{prefix}_Skani_report_final.csv')  # Save final result to CSV
    result_df.to_csv(result_file_path, index=False)


def summary(prefix, outdir, skani_genome_size):
    prefix = prefix.pop()
    outdir = outdir.pop()
    
    # Organize reports directory
    report_dir = str(outdir) + "/%s_Report" % prefix
    report_script_dir = str(outdir) + "/%s_Report/scripts" % prefix
    
    
    Coverage = pd.read_csv("results/%s/%s_Report/data/%s_Final_Coverage.txt" % (prefix, prefix, prefix), sep=',', header=0)
    Coverage.rename(columns = {'Sample_name':'Sample'}, inplace = True)

    #kraken = pd.read_csv("results/%s/%s_Report/data/%s_Kraken_report_final.csv" % (prefix, prefix, prefix), sep=',', header=0)
    
    mlst = pd.read_csv("results/%s/%s_Report/data/%s_MLST_results.csv" % (prefix, prefix, prefix), sep='\t', header=0)
    #mlst = mlst.replace(['_contigs_l1000.fasta'], '', regex=True)
    #mlst = mlst.replace(['results/.*/spades/'], '', regex=True)
    #mlst = mlst.replace(['%s' % prefix], '', regex=True)
    mlst['Sample'] = mlst['Sample'].replace(r'.*/spades/(.*?)/.*', r'\1', regex=True)

    multiqc_fastqc_summary = pd.read_csv("results/%s/%s_Report/multiqc/%s_QC_report_data/multiqc_fastqc.txt" % (prefix, prefix, prefix), sep='\t', header=0)
    patternDel = "_R2"
    filter = multiqc_fastqc_summary['Sample'].str.contains(patternDel)
    multiqc_fastqc_summary = multiqc_fastqc_summary[~filter]
    aftertrim_filter = multiqc_fastqc_summary['Sample'].str.contains("_R1_trim_paired")
    raw_multiqc_fastqc_summary = multiqc_fastqc_summary[~aftertrim_filter]
    raw_multiqc_fastqc_summary = raw_multiqc_fastqc_summary.replace(['_R1'], '', regex=True)
    
    aftertrim_multiqc_fastqc_summary = multiqc_fastqc_summary[aftertrim_filter]
    aftertrim_multiqc_fastqc_summary = aftertrim_multiqc_fastqc_summary.replace(['_R1_trim_paired'], '', regex=True)
    aftertrim_multiqc_fastqc_summary = aftertrim_multiqc_fastqc_summary.add_prefix('After_trim_')
    aftertrim_multiqc_fastqc_summary.rename(columns = {'After_trim_Sample':'Sample'}, inplace = True)

    multiqc_general_stats_summary = pd.read_csv("results/%s/%s_Report/multiqc/%s_QC_report_data/multiqc_general_stats.txt" % (prefix, prefix, prefix), sep='\t', header=0)
    quast_filter = multiqc_general_stats_summary['Sample'].str.contains("_contigs_l1000")
    multiqc_quast = multiqc_general_stats_summary[quast_filter]
    multiqc_quast = multiqc_quast.replace(['_contigs_l1000'], '', regex=True)
    
    if 'QUAST_mqc-generalstats-quast-N50' in multiqc_quast.columns and 'QUAST_mqc-generalstats-quast-Total_length' in multiqc_quast.columns:
        multiqc_quast = multiqc_quast[["Sample", "QUAST_mqc-generalstats-quast-N50", "QUAST_mqc-generalstats-quast-Total_length"]]
        multiqc_quast = multiqc_quast.rename(columns={"QUAST_mqc-generalstats-quast-N50": "N50", "QUAST_mqc-generalstats-quast-Total_length": "Total length"})
    elif 'N50' in multiqc_quast.columns and 'Total length' in multiqc_quast.columns:
        multiqc_quast = multiqc_quast[["Sample", "N50", "Total length"]]
    #multiqc_quast = multiqc_quast[["Sample", "N50", "Total length"]]

    contig_distribution = pd.read_csv("results/%s/%s_Report/multiqc/%s_QC_report_data/mqc_quast_num_contigs_1.txt" % (prefix, prefix, prefix), sep='\t', header=0)
    contig_distribution = contig_distribution.replace(['_contigs_l1000'], '', regex=True)
    contig_distribution['Total # of contigs'] = contig_distribution.sum(axis=1, numeric_only=True)
    contig_distribution = contig_distribution[['Sample', 'Total # of contigs']]
    
    #read final skani output file
    skani_summary = pd.read_csv("results/%s/%s_Report/data/%s_Skani_report_final.csv" % (prefix, prefix, prefix), sep=',', skipinitialspace=True, header=0, engine='python')

    QC_summary_temp1 = pd.merge(Coverage, mlst, on=["Sample", "Sample"],  how='left')
    QC_summary_temp2 = QC_summary_temp1
    QC_summary_temp3 = pd.merge(QC_summary_temp2, raw_multiqc_fastqc_summary, on=["Sample", "Sample"], how='left')
    QC_summary_temp4 = pd.merge(QC_summary_temp3, aftertrim_multiqc_fastqc_summary, on=["Sample", "Sample"], how='left')
    QC_summary_temp5 = pd.merge(QC_summary_temp4, multiqc_quast, on=["Sample", "Sample"], how='left')
    QC_summary_temp6 = pd.merge(QC_summary_temp5, contig_distribution, on=["Sample", "Sample"], how='left')
    
    #QC_summary_temp7 = QC_summary_temp6[["Sample" , "Total_reads" , "Total_bp" , "MeanReadLength" , "Coverage" , "Scheme" , "ST" , "PercentageofreadsforSpecies" , "#ofreadsforSpecies" , "Species" , "After_trim_per_base_sequence_content" , "After_trim_overrepresented_sequences" , "After_trim_%GC" , "After_trim_Total Bases" , "After_trim_Total Sequences" , "After_trim_median_sequence_length" , "After_trim_avg_sequence_length" , "After_trim_total_deduplicated_percentage" , "After_trim_Sequence length" , "After_trim_adapter_content" , "N50" , "Total length" , "Total # of contigs"]].copy() #.copy() to deal with SettingWithCopyWarning error
    QC_summary_temp7 = QC_summary_temp6[["Sample" , "Total_reads" , "Total_bp" , "MeanReadLength" , "Coverage" , "Scheme" , "ST" , "After_trim_per_base_sequence_content" , "After_trim_overrepresented_sequences" , "After_trim_%GC" , "After_trim_Total Bases" , "After_trim_Total Sequences" , "After_trim_median_sequence_length" , "After_trim_avg_sequence_length" , "After_trim_total_deduplicated_percentage" , "After_trim_Sequence length" , "After_trim_adapter_content" , "N50" , "Total length" , "Total # of contigs"]].copy() #.copy() to deal with SettingWithCopyWarning error

    skani_genome_size = list(skani_genome_size)[0]
    
    # Read in skani species gneome size table
    skani_genome_table = pd.read_csv(skani_genome_size)

    # Check the total length against the assembly length with ±15% rule
    def check_assembly_length(total_length, assembly_length):
        if pd.isnull(total_length):
            return 'FAIL'
        if pd.isnull(assembly_length):
            if config["genome_size"] <= total_length <= config["assembly_length"]:
                return 'PASS'
            else:
                return 'FAIL'
        lower_bound = assembly_length * 0.85
        upper_bound = assembly_length * 1.15
        if lower_bound <= total_length <= upper_bound:
            return 'PASS'
        else:
            return 'FAIL'
    
    QC_summary_temp8 = pd.merge(QC_summary_temp7, skani_summary, on=["Sample", "Sample"], how='left') # Merge skani df into the existing dataframe

    # Merge QC summary with skani_genome_size on species
    QC_summary_temp8 = QC_summary_temp8.merge(
        skani_genome_table,
        left_on='Species',  
        right_on='Species', 
        how='left'
    )

    # Apply the length check function
    QC_summary_temp8['Length Check'] = QC_summary_temp8.apply(
        lambda row: check_assembly_length(row['Total length'], row['Assembly_Length']),
        axis=1
    )

    # Updated QC check based on the new Length Check condition
    QC_check_condition = [
        (QC_summary_temp8['Total # of contigs'] > config["max_contigs"]),
        (QC_summary_temp8['Total # of contigs'] < config["min_contigs"]),
        (QC_summary_temp8['Coverage'] < config["coverage"]),
        (QC_summary_temp8['Length Check'] == 'FAIL'),
        (QC_summary_temp8['Total # of contigs'].isnull()),
        (pd.isnull(QC_summary_temp8['Total length'])),
    ]

    status = ['FAIL', 'FAIL', 'FAIL', 'FAIL', 'FAIL', "Run FAIL"]

    QC_summary_temp8['QC Check'] = np.select(QC_check_condition, status, default='PASS')
    
     # Remove the 'Length Check' column
    QC_summary_temp8 = QC_summary_temp8.drop(columns=['Length Check', 'Assembly_Length'])

    # First, get the current list of columns
    columns = list(QC_summary_temp8.columns)
    
    # Insert QC Check between Total # of contigs and ANIs
    # Need to know the index positions of these columns to make the rearrangement correctly
    contigs_index = columns.index('Total # of contigs')
    ani_index = columns.index('ANI')
    qc_check_index = columns.index('QC Check')

    # Create the new column order
    # Put all columns before Total # of contigs then Total # of contigs, QC Check, ANI and the rest
    new_columns = columns[:contigs_index + 1] + ['QC Check'] + columns[ani_index:qc_check_index]

    # Rearrange the columns
    #QC_summary_temp8 = QC_summary_temp8[new_columns]
    
    # Read in samples that failed coverage and assembly respectively
    failed_coverage_path = f"results/{prefix}/sample_files/samples_failed_coverage_summary.csv"
    failed_assembly_path = f"results/{prefix}/sample_files/samples_failed_assembly_summary.csv"

    # Initialize empty DataFrames for optional files
    failed_coverage_df = pd.DataFrame()
    failed_assembly_df = pd.DataFrame()

    # Read existing files
    if os.path.exists(failed_coverage_path):
        failed_coverage_df = pd.read_csv(failed_coverage_path)

    if os.path.exists(failed_assembly_path):
        failed_assembly_df = pd.read_csv(failed_assembly_path)

    # Find common columns with QC_summary_temp8
    if not failed_coverage_df.empty or not failed_assembly_df.empty:
        common_columns = QC_summary_temp8.columns
        if not failed_coverage_df.empty:
            common_columns = common_columns.intersection(failed_coverage_df.columns)
        if not failed_assembly_df.empty:
            common_columns = common_columns.intersection(failed_assembly_df.columns)

        # Select only the matching columns
        QC_summary_temp8 = QC_summary_temp8[common_columns]
        if not failed_coverage_df.empty:
            failed_coverage_df = failed_coverage_df[common_columns]
        if not failed_assembly_df.empty:
            failed_assembly_df = failed_assembly_df[common_columns]

        # Concatenate available DataFrames
        QC_summary_temp8 = pd.concat([QC_summary_temp8, failed_coverage_df, failed_assembly_df], ignore_index=True)

        # Fill missing values
        QC_summary_temp8 = QC_summary_temp8.fillna("NA")
        QC_summary_temp8 = QC_summary_temp8[new_columns]
        
    else:
        # Rearrange the columns if no failed coverage/assembly files exist
        QC_summary_temp8 = QC_summary_temp8[new_columns]

        # # Ensure all dataframes have the same columns as QC_summary_temp8
        # common_columns = QC_summary_temp8.columns.intersection(failed_assembly_df.columns).intersection(failed_coverage_df.columns)
        
        # # Select only the matching columns
        # QC_summary_temp8 = QC_summary_temp8[common_columns]
        # failed_assembly_df = failed_assembly_df[common_columns]
        # failed_coverage_df = failed_coverage_df[common_columns]

        # # Concatenate the dataframes
        # QC_summary_temp8 = pd.concat([QC_summary_temp8, failed_assembly_df, failed_coverage_df], ignore_index=True)
        
        # QC_summary_temp8 = QC_summary_temp8.fillna("NA")
    
    # Remove duplicates based on the "Sample" column
    if "Sample" in QC_summary_temp8.columns:
        QC_summary_temp8 = QC_summary_temp8.drop_duplicates(subset=["Sample"])

    QC_summary_temp8.to_csv('results/%s/%s_Report/data/%s_QC_summary.csv' % (prefix, prefix, prefix), index=False)

def plot(prefix, outdir):
    prefix = prefix.pop()
    outdir = outdir.pop()
    
    # Organize reports directory
    report_dir = str(outdir) + "/%s_Report" % prefix
    report_script_dir = str(outdir) + "/%s_Report/scripts" % prefix
    
    QC_summary = pd.read_csv('results/%s/%s_Report/data/%s_QC_summary.csv' % (prefix, prefix, prefix), sep=',', header=0)    

    Coverage = pd.read_csv("results/%s/%s_Report/data/%s_Final_Coverage.txt" % (prefix, prefix, prefix), sep=',', header=0)    
    Coverage_dist = QC_summary.sort_values(by='Coverage',ascending=False).plot(x='Sample', y='Coverage', kind="barh", title="Estimated Genome Coverage", figsize=(20, 20), fontsize=40).get_figure()
    Coverage_dist.savefig('%s/fig/%s_Coverage_distribution.png' % (report_dir, prefix), dpi=600)


    ax1 = QC_summary.plot.scatter(x = 'After_trim_total_deduplicated_percentage', y = 'After_trim_Total Sequences', c = 'DarkBlue')
    fig = ax1.get_figure()
    fig.savefig('%s/fig/%s_raw_dedup_vs_totalsequence.png' % (report_dir, prefix), dpi=600)

    ax1 = QC_summary.plot.scatter(x = 'After_trim_total_deduplicated_percentage', y = 'After_trim_Total Sequences', c = 'DarkBlue')
    fig = ax1.get_figure()
    fig.savefig('%s/fig/%s_aftertrim_dedup_vs_totalsequence.png' % (report_dir, prefix), dpi=600)
    ax1.cla()


rule all:
    input:
        coverage_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Final_Coverage.txt", prefix=PREFIX),
        #kraken_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Kraken_report_final.csv", prefix=PREFIX),
        skani_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Skani_report_final.csv", prefix=PREFIX),
        multiqc_report = expand("results/{prefix}/{prefix}_Report/multiqc/{prefix}_QC_report.html", prefix=PREFIX),
        mlst_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_MLST_results.csv", prefix=PREFIX),
        QC_summary = expand("results/{prefix}/{prefix}_Report/data/{prefix}_QC_summary.csv", prefix=PREFIX),
        #QC_plot = expand("results/{prefix}/{prefix}_Report/fig/{prefix}_Coverage_distribution.png", prefix=PREFIX)

rule coverage_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        coverage_out = expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", prefix=PREFIX, sample=SAMPLE_PASS_ASSEMBLY)
    output:
        coverage = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_Final_Coverage.txt",
    params:
        prefix = "{prefix}",
    run:
        coverage_report({params.prefix}, {input.outdir})

rule skani_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        skani_out = expand("results/{prefix}/skani/{sample}/{sample}_skani_output.txt", prefix=PREFIX, sample=SAMPLE_PASS_ASSEMBLY)
    output:
        skani_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_Skani_report_final.csv",
    params:
        prefix = "{prefix}",
    run:
        skani_report({input.outdir}, {params.prefix})

rule multiqc:
    input:
        inputdir = lambda wildcards: expand(f"results/{wildcards.prefix}"),
        coverage = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Final_Coverage.txt"),
        #kraken = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Kraken_report_final.csv"),
        mlst = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_MLST_results.csv"),
    output:
        multiqc_fastqc_report = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report.html",
        multiqc_fastqc = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report_data/multiqc_fastqc.txt",
        multiqc_general_stats = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report_data/multiqc_general_stats.txt",
        #fastqc_report = f"results//{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report_data/multiqc_fastqc.txt"
    params:
        outdir = "results/{prefix}/{prefix}_Report",
        prefix = "{prefix}",
    #conda:
    #    "envs/multiqc.yaml"
    singularity:
        "docker://staphb/multiqc:1.19"
    shell:
        "multiqc -f --export --outdir {params.outdir}/multiqc -n {params.prefix}_QC_report -i {params.prefix}_QC_report {input.inputdir}/quality_aftertrim/*/*_Forward {input.inputdir}/prokka/* {input.inputdir}/quast/*"

rule mlst_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        mlst_out = expand("results/{prefix}/mlst/{sample}/report.tsv", prefix=PREFIX, sample=SAMPLE_PASS_ASSEMBLY)
    output:
        mlst_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_MLST_results.csv",
    params:
        prefix = "{prefix}",
    shell:
        "echo \"Sample\tScheme\tST\" > {output.mlst_report} && cut -f1-3 {input.outdir}/mlst/*/report.tsv >> {output.mlst_report}"

rule Summary:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        multiqc_fastqc_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/multiqc/{wildcards.prefix}_QC_report.html"),
        multiqc_fastqc = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/multiqc/{wildcards.prefix}_QC_report_data/multiqc_fastqc.txt"),
        multiqc_general_stats = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/multiqc/{wildcards.prefix}_QC_report_data/multiqc_general_stats.txt"),
        coverage = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Final_Coverage.txt"),
        #kraken = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Kraken_report_final.csv"),
        mlst = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_MLST_results.csv"),
        skani_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Skani_report_final.csv"),
    output:
        QC_summary_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_QC_summary.csv",
    params:
        prefix = "{prefix}",
        skani_genome_size_table = config["skani_genome_size"] # In the config folder
    run:
        summary({params.prefix}, {input.outdir}, {params.skani_genome_size_table})

rule plot:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        QC_summary_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_QC_summary.csv"),
    output:
        QC_summary_report = f"results/{{prefix}}/{{prefix}}_Report/fig/{{prefix}}_Coverage_distribution.png",
    params:
        prefix = "{prefix}",
    run:
        plot({params.prefix}, {input.outdir})
