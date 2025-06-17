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
    
rule skani_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        skani_out = lambda wildcards: [
            f"results/{prefix}/skani/{sample}/{sample}_skani_output.txt"
            for sample in samples_that_passed_assembly(wildcards, return_samples_only=True)
        ] # expand("results/{prefix}/skani/{sample}/{sample}_skani_output.txt", prefix=PREFIX, sample=SAMPLE)
    output:
        skani_report = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_Skani_report_final.csv",
    params:
        prefix = "{prefix}",
    run:
        skani_report({input.outdir}, {params.prefix})
