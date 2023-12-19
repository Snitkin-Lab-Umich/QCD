# Author: Ali Pirani
configfile: "config/config.yaml"

import pandas as pd
import os
import json
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

PREFIX = config["prefix"]

SHORTREADS = list(samples_df['sample_id'])

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)

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

    for sampl in SAMPLE:
        coverage_json = "results/%s/%s/raw_coverage/%s_coverage.json" % (prefix, sampl, sampl)
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

def kraken_report(prefix, outdir):
    prefix = prefix.pop()
    outdir = outdir.pop()
    
    # Organize reports directory
    report_dir = str(outdir) + "/%s_Report" % prefix
    report_script_dir = str(outdir) + "/%s_Report/scripts" % prefix

    kraken_dir = str(outdir) + "/*/kraken"

    kraken_summary_script = open("%s/kraken_summary.sh" % report_script_dir, 'w+')
    kraken_summary_script.write("echo \"Sample,Percentage of reads for Species,# of reads for Species, Species\" > %s/data/%s_Kraken_report_final.csv\n" % (report_dir, prefix))
    kraken_summary_script.write("for i in results/%s/*/kraken/*_kraken_report.tsv; do grep -w 'S' $i | sort -k1n | tail -n1; done > /tmp/Kraken_report_temp.txt\n" % prefix)
    kraken_summary_script.write("ls -d results/%s/*/kraken/*_kraken_report.tsv | awk -F'/' '{print $NF}' | sed 's/_kraken_report.tsv//g' > %s/data/samplenames.txt\n" % (prefix, report_dir))
    kraken_summary_script.write("paste %s/data/samplenames.txt /tmp/Kraken_report_temp.txt > /tmp/Kraken_report_combined.txt\n" % (report_dir))
    kraken_summary_script.write("awk -F'\\t' 'BEGIN{OFS=\",\"};{print $1, $2, $3, $7}' /tmp/Kraken_report_combined.txt >> %s/data/%s_Kraken_report_final.csv\n" % (report_dir, prefix))
    kraken_summary_script.write("sed -i 's/\s//g' %s/data/%s_Kraken_report_final.csv\n" % (report_dir, prefix))
    kraken_summary_script.close()

    os.system("bash %s/kraken_summary.sh" % report_script_dir)

def summary(prefix, outdir):
    prefix = prefix.pop()
    outdir = outdir.pop()
    
    # Organize reports directory
    report_dir = str(outdir) + "/%s_Report" % prefix
    report_script_dir = str(outdir) + "/%s_Report/scripts" % prefix
    
    
    Coverage = pd.read_csv("results/%s/%s_Report/data/%s_Final_Coverage.txt" % (prefix, prefix, prefix), sep=',', header=0)
    Coverage.rename(columns = {'Sample_name':'Sample'}, inplace = True)

    kraken = pd.read_csv("results/%s/%s_Report/data/%s_Kraken_report_final.csv" % (prefix, prefix, prefix), sep=',', header=0)
    
    mlst = pd.read_csv("results/%s/%s_Report/data/%s_MLST_results.csv" % (prefix, prefix, prefix), sep='\t', header=0)
    mlst = mlst.replace(['_contigs_l1000.fasta'], '', regex=True)
    mlst = mlst.replace(['results/.*/spades/'], '', regex=True)
    mlst = mlst.replace(['%s' % prefix], '', regex=True)

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

    QC_summary_temp1 = pd.merge(Coverage, mlst, on=["Sample", "Sample"],  how='left')
    QC_summary_temp2 = pd.merge(QC_summary_temp1, kraken, on=["Sample", "Sample"],  how='left')
    QC_summary_temp3 = pd.merge(QC_summary_temp2, raw_multiqc_fastqc_summary, on=["Sample", "Sample"], how='left')
    QC_summary_temp4 = pd.merge(QC_summary_temp3, aftertrim_multiqc_fastqc_summary, on=["Sample", "Sample"], how='left')
    QC_summary_temp5 = pd.merge(QC_summary_temp4, multiqc_quast, on=["Sample", "Sample"], how='left')
    QC_summary_temp6 = pd.merge(QC_summary_temp5, contig_distribution, on=["Sample", "Sample"], how='left')
    
    QC_summary_temp7 = QC_summary_temp6[["Sample" , "Total_reads" , "Total_bp" , "MeanReadLength" , "Coverage" , "Scheme" , "ST" , "PercentageofreadsforSpecies" , "#ofreadsforSpecies" , "Species" , "After_trim_per_base_sequence_content" , "After_trim_overrepresented_sequences" , "After_trim_%GC" , "After_trim_Total Bases" , "After_trim_Total Sequences" , "After_trim_median_sequence_length" , "After_trim_avg_sequence_length" , "After_trim_total_deduplicated_percentage" , "After_trim_Sequence length" , "After_trim_adapter_content" , "N50" , "Total length" , "Total # of contigs"]]

    QC_check_condition = [
    (QC_summary_temp7['Total # of contigs'] > config["max_contigs"]),
    (QC_summary_temp7['Total # of contigs'] < config["min_contigs"]),
    (QC_summary_temp7['Total length'] > config["assembly_length"]),
    (QC_summary_temp7['Coverage'] < config["coverage"]),
    (QC_summary_temp7['Total # of contigs'].isnull()),
    ]

    status = ['FAIL', 'FAIL', 'FAIL', 'FAIL', "Run FAIL"]

    QC_summary_temp7['QC Check'] = np.select(QC_check_condition, status)

    QC_summary_temp7.to_csv('results/%s/%s_Report/data/%s_QC_summary.csv' % (prefix, prefix, prefix), index=False)

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

    ax = sns.scatterplot(x=QC_summary['Total # of contigs'], y=QC_summary['After_trim_%GC'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_Aftertrim_GC.png' % (report_dir, prefix), dpi=600)
    plt.savefig('%s/fig/%s_Assembly_contig_vs_Aftertrim_GC.png' % (report_dir, prefix), dpi=200)
    ax.cla()

    ax = sns.scatterplot(x=QC_summary['Total length'], y=QC_summary['After_trim_%GC'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_Aftertrim_GC.png' % (report_dir, prefix), dpi=600)
    plt.savefig('%s/fig/%s_Assembly_length_vs_Aftertrim_GC.png' % (report_dir, prefix), dpi=200)
    ax.cla()

    ax = sns.scatterplot(x=QC_summary['Total # of contigs'], y=QC_summary['N50'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_N50.png' % (report_dir, prefix), dpi=600)
    plt.savefig('%s/fig/%s_Assembly_contig_vs_N50.png' % (report_dir, prefix), dpi=200)
    ax.cla()

    ax = sns.scatterplot(x=QC_summary['Total # of contigs'], y=QC_summary['Coverage'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_N50.png' % (report_dir, prefix), dpi=600)
    plt.savefig('%s/fig/%s_Assembly_contig_vs_Coverage.png' % (report_dir, prefix), dpi=200)
    ax.cla()

    ax = sns.scatterplot(x=QC_summary['Total # of contigs'], y=QC_summary['Total length'], hue=QC_summary['Species'], s=100, style=QC_summary['Species'])
    #g.legend(loc='right', bbox_to_anchor=(1.30, 0.5), ncol=1)
    #fig2 = g.get_figure()
    #fig2.savefig('%s/fig/%s_Assembly_contig_vs_N50.png' % (report_dir, prefix), dpi=600)
    plt.savefig('%s/fig/%s_Assembly_contig_vs_length.png' % (report_dir, prefix), dpi=200)
    ax.cla()


rule all:
    input:
        coverage_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Final_Coverage.txt", prefix=PREFIX),
        kraken_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_Kraken_report_final.csv", prefix=PREFIX),
        mlst_report = expand("results/{prefix}/{prefix}_Report/data/{prefix}_MLST_results.csv", prefix=PREFIX),
        multiqc_report = expand("results/{prefix}/{prefix}_Report/multiqc/{prefix}_QC_report.html", prefix=PREFIX),
        QC_summary = expand("results/{prefix}/{prefix}_Report/data/{prefix}_QC_summary.csv", prefix=PREFIX),
        QC_plot = expand("results/{prefix}/{prefix}_Report/fig/{prefix}_Coverage_distribution.png", prefix=PREFIX),

rule coverage_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
    output:
        coverage = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_Final_Coverage.txt",
    params:
        prefix = "{prefix}",
    run:
        coverage_report({params.prefix}, {input.outdir})

rule amr_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
    output:
        amr_summary = f"results/{{prefix}}/report/{{prefix}}_AMR_minimal_report.csv",
    params:
        prefix = "{prefix}",
        phandango = "--no_tree"
    conda:
        "envs/ariba.yaml"
    shell:
        "ariba summary --preset minimal {params.phandango} {input.outdir}/report/{params.prefix}_AMR_minimal_report {input.outdir}/*/ariba_card/report.tsv && ariba summary --preset all {params.phandango} {input.outdir}/report/{params.prefix}_AMR_all_report {input.outdir}/*/ariba_card/report.tsv"

rule kraken_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
    output:
        kraken_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_Kraken_report_final.csv",
    params:
        prefix = "{prefix}",
    run:
        kraken_report({params.prefix}, {input.outdir})

rule multiqc:
    input:
        inputdir = lambda wildcards: expand(f"results/{wildcards.prefix}"),
        coverage = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Final_Coverage.txt"),
        kraken = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Kraken_report_final.csv"),
        mlst = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_MLST_results.csv"),
    output:
        multiqc_fastqc_report = f"results/{{prefix}}/{{prefix}}_Report/multiqc/{{prefix}}_QC_report.html",
    params:
        outdir = "results/{prefix}/{prefix}_Report",
        prefix = "{prefix}",
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc -f --export --outdir {params.outdir}/multiqc -n {params.prefix}_QC_report -i {params.prefix}_QC_report {input.inputdir}/*/quality_aftertrim/*_Forward {input.inputdir}/*/kraken {input.inputdir}/*/prokka {input.inputdir}/*/quast"

rule mlst:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
    output:
        mlst_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_MLST_results.csv",
    params:
        prefix = "{prefix}",
    shell:
        "echo \"Sample\tScheme\tST\" > {output.mlst_report} && cut -f1-3 {input.outdir}/*/mlst/report.tsv >> {output.mlst_report}"

rule Summary:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        multiqc_fastqc_report = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/multiqc/{wildcards.prefix}_QC_report.html"),
        coverage = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Final_Coverage.txt"),
        kraken = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_Kraken_report_final.csv"),
        mlst = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.prefix}_Report/data/{wildcards.prefix}_MLST_results.csv"),
    output:
        QC_summary_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_QC_summary.csv",
    params:
        prefix = "{prefix}",
    run:
        summary({params.prefix}, {input.outdir})

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
