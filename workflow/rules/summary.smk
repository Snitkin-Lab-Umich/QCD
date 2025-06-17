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
    QC_summary_temp8 = QC_summary_temp8[new_columns]

    QC_summary_temp8 = QC_summary_temp8.fillna("NA")
    
    QC_summary_temp8.to_csv('results/%s/%s_Report/data/%s_QC_summary.csv' % (prefix, prefix, prefix), index=False)
    
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