
def coverage_report(prefix, outdir):
    prefix = prefix.pop()
    report_dir = str(outdir.pop()) + "/%s_Report" % prefix
    # Generate Coverage report 
    final_coverage_file = "%s/data/%s_Final_Coverage.txt" % (report_dir, prefix)
    f3=open(final_coverage_file, 'w+')
    header = "Sample,Total_reads,Total_bp,MeanReadLength,Coverage\n"
    f3.write(header)

    for sampl in SAMPLE:
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

rule coverage_report:
    input:
        outdir = lambda wildcards: expand(f"results/{wildcards.prefix}/"),
        coverage_out = expand("results/{prefix}/raw_coverage/{sample}/{sample}_coverage.json", prefix=PREFIX, sample=SAMPLE)
    output:
        coverage = f"results/{{prefix}}/{{prefix}}_Report/data/{{prefix}}_Final_Coverage.txt",
    params:
        prefix = "{prefix}",
    run:
        coverage_report({params.prefix}, {input.outdir})