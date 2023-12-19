# Import Libraries and set path
import numpy as np
import pandas as pd
from IPython.display import HTML
import os
import readline
import argparse
from itertools import islice
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import argparse
import json

def parser():
    parser = argparse.ArgumentParser(description='\nGather outputs and summarize QC reports\n')
    parser.add_argument('-prefix', action='store', dest="prefix", help='Analysis prefix name for QC report', required=False)
    parser.add_argument('-outdir', action='store', dest="outdir", help='Path to QCD results folder.', required=True)
    parser.add_argument('-samples', action='store', dest="samples", help='List of sample names', required=True, default="PE")
    return parser
args = parser().parse_args()

report_dir = args.outdir + "/report"
report_script_dir = args.outdir + "/report/scripts"
isExist = os.path.exists(report_dir)
if not isExist:
   os.makedirs(report_dir)
isExist = os.path.exists(report_script_dir)
if not isExist:
   os.makedirs(report_script_dir)
   
# Generate Coverage report 
print("- Generating Coverage report.")
final_coverage_file = "%s/%s_Final_Coverage.txt" % (report_dir, args.prefix)
f3=open(final_coverage_file, 'w+')
header = "Sample_name,Total_reads,Total_bp,MeanReadLength,Coverage\n"
f3.write(header)
with open(args.samples, 'r') as samples:
    for line in samples:
        line = line.strip()
        sampl = line.split(',')
        
        if not sampl[0].startswith('sample_id'):
            coverage_json = "%s/%s/raw_coverage/%s_coverage.json" % (args.outdir, sampl[0], sampl[0])
            f = open(coverage_json)
            data = json.load(f)
            # data = json.loads(coverage_json)
            f3.write("%s,%s,%s,%s,%s\n" % (sampl[0], data['qc_stats']['read_total'], data['qc_stats']['total_bp'], data['qc_stats']['read_mean'], data['qc_stats']['coverage']))

f3.close()
Coverage = pd.read_csv(final_coverage_file, sep=',', header=0)
Coverage = Coverage.replace(['_R1.fastq.gz'], '', regex=True)

print ("Number of Samples in Coverage Report - %s" % len(Coverage))

Coverage_NEG_CNTL = Coverage[Coverage.Sample_name.str.match('(.*NEG*)')]

print ("Number of Negative Control samples %s" % len(Coverage_NEG_CNTL))

print ("Number of Negative Control samples with > 100X coverage %s" % len(Coverage_NEG_CNTL[Coverage_NEG_CNTL['Coverage'] > 100]))

Coverage_dist = Coverage.sort_values(by='Coverage',ascending=False).plot(x='Sample_name', y='Coverage', kind="barh", title="Estimated Genome Coverage", figsize=(20, 20), fontsize=40).get_figure()

Coverage_dist.savefig('%s/%s_Coverage_distribution.pdf' % (report_dir, args.prefix))

# print("- Generating Ariba AMR report.")
# ariba_dir = args.outdir + "/*/ariba_card/"

# ariba_minimal_summary = "/nfs/turbo/umms-esnitkin/conda/ariba/bin/ariba summary --preset minimal %s_AMR_minimal_report %s/report.tsv" % (args.prefix, ariba_dir)
# ariba_all_summary = "/nfs/turbo/umms-esnitkin/conda/ariba/bin/ariba summary --preset minimal %s_AMR_all_report %s/report.tsv" % (args.prefix, ariba_dir)
# ariba_summary_script = open("%s/ariba_summary.sh" % report_script_dir, 'w+')
# ariba_summary_script.write(ariba_minimal_summary + '\n' + ariba_all_summary)
# ariba_summary_script.close()

# os.system(ariba_minimal_summary)
# os.system(ariba_all_summary)

print("- Generating Kraken report.")
kraken_dir = args.outdir + "/*/kraken"

kraken_summary_script = open("%s/kraken_summary.sh" % report_script_dir, 'w+')
kraken_summary_script.write("echo \"Sample,Percentage of reads for Species,# of reads for Species, Species\" > %s/Kraken_report_final.csv\n" % report_dir)
kraken_summary_script.write("for i in %s/*_kraken_report.tsv; do grep -w 'S' $i | sort -k1n | tail -n1; done > /tmp/Kraken_report_temp.txt\n" % (kraken_dir))
kraken_summary_script.write("ls -d %s/*_kraken_report.tsv | awk -F'/' '{print $NF}' | sed 's/_kraken_report.tsv//g' > %s/samplenames.txt\n" % (kraken_dir, report_dir))
kraken_summary_script.write("paste %s/samplenames.txt /tmp/Kraken_report_temp.txt > /tmp/Kraken_report_combined.txt\n" % (report_dir))
kraken_summary_script.write("awk -F'\\t' 'BEGIN{OFS=\",\"};{print $1, $2, $3, $7}' /tmp/Kraken_report_combined.txt >> %s/Kraken_report_final.csv\n" % (report_dir))
kraken_summary_script.write("sed -i 's/\s//g' %s/Kraken_report_final.csv\n" % report_dir)
kraken_summary_script.close()

os.system("%s/kraken_summary.sh" % report_script_dir)

