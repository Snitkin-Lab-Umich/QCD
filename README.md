# QCD - Quality Control and Contamination Detection workflow.

QCD is a workflow for microbial Illumina sequencing quality control and contamination detection implemented in Snakemake.

### Summary

As part of the SOP in the [Snitkin lab](https://thesnitkinlab.com/index.php), this pipeline should be run on raw sequencing data as soon as the data is available from the sequencing core. 

In short, it performs the following steps:

* [Fastqc](https://github.com/s-andrews/FastQC) is used to generate HTML reports to asses quality of sequencing reads before and after trimming reads. 
* Trims and filters low-quality bases and adapter sequences from raw FASTQ reads using [Trimmomatic](https://github.com/usadellab/Trimmomatic).
* [fastq-scan](https://github.com/rpetit3/fastq-scan) is used to estimate genome coverage of FASTQ files.
* Assembles trimmed reads into contigs using [SPAdes](https://github.com/ablab/spades).
* The assembled contigs from [SPAdes](https://github.com/ablab/spades) is then passed through [Prokka](https://github.com/tseemann/prokka) for annotation, [QUAST](https://quast.sourceforge.net/) for assembly statistics, [MLST](https://github.com/tseemann/mlst) for determining sequence type based on sequences of housekeeping genes, [skani](https://github.com/bluenote-1577/skani) to identify closest reference genome and [BUSCO](https://busco.ezlab.org/) for assembly completeness statistics.
* [Multiqc](https://github.com/MultiQC/MultiQC) aggregates the final outputs from [Fastqc](https://github.com/s-andrews/FastQC), [Prokka](https://github.com/tseemann/prokka) and [QUAST](https://quast.sourceforge.net/) to produce a HTML report.
* The final step in the pipeline is to generate a QC report that will specify which samples have passed or failed the pipeline.

The workflow generates all the output in the output prefix folder set in the config file (instructions on setup found [below](#config)). Each workflow steps gets its own individual folder as shown. **Note that this overview does not capture all possible outputs from each tool; it only highlights the primary directories and some of their contents.**

```
results/2025-05-29_Project_MRSA_QCD/
├── 2025-05-29_Project_MRSA_QCD_Report
│   ├── data
│   │   ├── 2025-05-29_Project_MRSA_QCD_Final_Coverage.txt
│   │   ├── 2025-05-29_Project_MRSA_QCD_MLST_results.csv
│   │   ├── 2025-05-29_Project_MRSA_QCD_QC_summary.csv
│   │   └── 2025-05-29_Project_MRSA_QCD_Skani_report_final.csv
│   ├── fig
│   ├── multiqc
├── downsample
│   ├── MRSA_jail_100
│   │   ├── MRSA_jail_100_R1_trim_paired.fastq.gz
│   │   └── MRSA_jail_100_R2_trim_paired.fastq.gz
│   └── MRSA_jail_275
│       ├── MRSA_jail_275_R1_trim_paired.fastq.gz
│       └── MRSA_jail_275_R2_trim_paired.fastq.gz
├── mlst
│   └── MRSA_jail_100
│       └── report.tsv
├── prokka
│   └── MRSA_jail_100
│       ├── MRSA_jail_100.gbk
│       ├── MRSA_jail_100.gff
│       └── MRSA_jail_100.txt
├── quality_aftertrim
│   ├── MRSA_jail_100
│   │   ├── MRSA_jail_100_Forward
│   │   │   ├── MRSA_jail_100_R1_trim_paired_fastqc.html
│   │   │   └── MRSA_jail_100_R1_trim_paired_fastqc.zip
├── quality_raw
│   ├── MRSA_jail_100
│   │   ├── MRSA_jail_100_Forward
│   │   │   ├── MRSA_jail_100_R1_fastqc.html
│   │   │   └── MRSA_jail_100_R1_fastqc.zip
├── quast
│   └── MRSA_jail_100
│       └── transposed_report.txt
├── raw_coverage
│   ├── MRSA_jail_100
│   │   └── MRSA_jail_100_coverage.json
├── skani
│   └── MRSA_jail_100
│       └── MRSA_jail_100_skani_output.txt
├── spades
│   ├── MRSA_jail_100
│   │   ├── contigs.fasta
│   │   └── MRSA_jail_100_contigs_l1000.fasta
└── trimmomatic
    ├── MRSA_jail_100
    │   ├── MRSA_jail_100_R1_trim_paired.fastq.gz
    │   └── MRSA_jail_100_R2_trim_paired.fastq.gz

```


## Installation 


> If you are using Great Lakes HPC, ensure you are cloning the repository in your scratch directory. Change `your_uniqname` to your uniqname. 

```

cd /scratch/esnitkin_root/esnitkin1/your_uniqname/

```

> Clone the github directory onto your system. 

```

git clone https://github.com/Snitkin-Lab-Umich/QCD.git

```

> Ensure you have successfully cloned QCD. Type `ls` and you should see the newly created directory **_QCD_**. Move to the newly created directory.

```

cd QCD

```

> Load Bioinformatics, snakemake and singularity modules from Great Lakes modules.

```

module load Bioinformatics snakemake singularity

```
<!--
```

module load snakemake singularity

```
-->

This workflow makes use of singularity containers available through [State Public Health Bioinformatics group](https://github.com/StaPH-B/docker-builds). If you are working on Great Lakes (umich cluster)—you can load snakemake and singularity modules as shown above. However, if you are running it on your local or other computing platform, ensure you have snakemake and singularity installed.


## Setup config, samples and cluster files

**_If you are just testing this pipeline, the config and sample files are already loaded with test data, so you do not need to make any additional changes to them. However, it is a good idea to change the prefix (name of your output folder) in the config file to give you an idea of what variables need to be modified when running your own samples on QCD._**

### Config
As an input, the snakemake file takes a config file where you can set the path to `samples.csv`, path to your raw sequencing reads, path to adapter fasta file etc. Instructions on how to modify `config/config.yaml` is found in `config.yaml`. 

### Samples
Add samples to `config/samples.csv` following the explanation provided below. `samples.csv` should be a comma seperated file consisting of two columns—`sample_id` and `illumina_r1`.

* `sample_id` is the prefix that should be extracted from your FASTQ reads. For example, in  your raw FASTQ files directory, if you have a file called `Rush_KPC_110_R1.fastq.gz`, your sample_id would be `Rush_KPC_110`.

* `illumina_r1` is the name of the entire raw FASTQ file. In the same directory,  if your file is called `Rush_KPC_110_R1.fastq.gz`, your sample_id would be `Rush_KPC_110_R1.fastq.gz`. **_Only include forward reads._**

You can create samples.csv file using the following for loop. Replace *path_to_your_raw_reads* below with the actual path to your raw sequencing reads.

```

echo "sample_id,illumina_r1" > config/samples.csv

for read1 in path_to_your_raw_reads/*_R1.fastq.gz; do
    sample_id=$(basename $read1 | sed 's/_R1.fastq.gz//g')
    read1_basename=$(basename $read1)
    echo $sample_id,$read1_basename
done >> config/samples.csv

```

### Cluster file

Reduce the walltime (to ~6 hours) in `config/cluster.json` to ensure the jobs are being submitted in a timely manner. 

## Quick start

>  Start an interactive session before you run the pipeline. 

```
salloc --mem-per-cpu=10G --account=esnitkin1
```

### Run QCD on a set of samples.

> Preview the steps in QCD by performing a dryrun of the pipeline. 

```

snakemake -s workflow/QCD.smk --dryrun -p

```

> Run QCD locally

```

snakemake -s workflow/QCD.smk -p --configfile config/config.yaml --cores all

```

>Run QCD directly on terminal (**_note: if you close your computer/shut down terminal, the pipeline will stop running. Terminal window has to be open until pipeline runs to completion._**)

```

snakemake -s workflow/QCD.smk -p --use-conda --use-singularity --use-envmodules -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock

```
> Submit QCD as a batch job (**reccommended**)

Change these `SBATCH` commands: `--job-name` to a more descriptive name like run_QCD, `--mail-user` to your email address, `--time` depending on the number of samples you have (should be more than what you specified in `cluster.json`). Feel free to make changes to the other flags if you are comfortable doing so. Once you have made the necessary changes, save the below script as `run_QCD.sbat`. Don't forget to submit QCD to Slurm! `sbatch bash_script_to_run_QCD.sbat`.

```
#!/bin/bash

#SBATCH --job-name=QCD
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=youremail@umich.edu
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=12:00:00
#SBATCH --account=esnitkin1
#SBATCH --partition=standard

# Load necessary modules
module load Bioinformatics
module load snakemake singularity

# Extract prefix from the YAML config file
PREFIX=$(grep '^prefix:' config/config.yaml | awk '{print $2}')

# Define the file paths dynamically
PASS_COV_FILE="results/${PREFIX}/sample_files/samples_passed_coverage.csv"
PASS_ASSEMBLY_FILE="results/${PREFIX}/sample_files/samples_passed_assembly.csv"

# Run Snakemake the first time -- until coverage
snakemake -s workflow/QCD.smk -p --use-conda --use-singularity --use-envmodules -j 999 \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" \
    --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock 

# samples_passed_coverage.csv should have been created in the Snakemake command above
# If not found, throw an error and exit
if [ ! -s "$PASS_COV_FILE" ] || [ "$(wc -l < "$PASS_COV_FILE")" -le 1 ]; then
    echo "Error: $PASS_COV_FILE is missing or does not have any samples. Exiting."
    exit 1
else
    echo "$PASS_COV_FILE detected. Running second part of the workflow."
fi

# Run Snakemake again to run the second part of the workflow --until assembly
snakemake -s workflow/QCD.smk -p --use-conda --use-singularity --use-envmodules -j 999 \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" \
    --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock 

# samples_passed_assembly.csv should have been created in the Snakemake command above
# If not found, throw an error and exit
if [ ! -s "$PASS_ASSEMBLY_FILE" ] || [ "$(wc -l < "$PASS_ASSEMBLY_FILE")" -le 1 ]; then
    echo "Error: $PASS_ASSEMBLY_FILE is missing or does not have any samples. Exiting."
    exit 1
else
    echo "$PASS_ASSEMBLY_FILE detected. Running second part of the workflow."
fi

# Run Snakemake to finish running the rest of the pipeline
snakemake -s workflow/QCD.smk -p --use-conda --use-singularity --use-envmodules -j 999 \
    --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" \
    --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 1000 --nolock 

# Run Snakemake for the last time to generate QC report 
snakemake -s workflow/QCD_report.smk -p --use-singularity --cores all

```

![Alt text](images/QCD_dag.svg)
<!--
### Gather Summary files and generate a report. 

>Start an interactive session in your current directory i.e. `QCD`.

```

srun --account=esnitkin1 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=5GB --cpus-per-task=1 --time=12:00:00 --pty /bin/bash

```

> Preview the steps in QCD report by performing a dryrun of the pipeline. 

```

snakemake -s QCD_report.smk --dryrun -p

```
> Run QCD report on Great lakes HPC

```

snakemake -s QCD_report.smk -p --use-singularity --cores 2

```
![Alt text](images/QCD_dag.svg)
-->
## Dependencies

### Near Essential
* [Snakemake>=7.32.4](https://snakemake.readthedocs.io/en/stable/#)
* [Conda](https://docs.conda.io/en/latest/)

<!--All the necessary software stack required for the workflow will be installed using conda package manager.-->

### Tool stack used in workflow

* [fastq-scan](https://github.com/rpetit3/fastq-scan)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [SPades](https://github.com/ablab/spades)
* [skani](https://github.com/bluenote-1577/skani)
* [bioawk](https://github.com/lh3/bioawk)
* [Prokka](https://github.com/tseemann/prokka)
* [mlst](https://github.com/tseemann/mlst)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [BUSCO](https://busco.ezlab.org/)
* [Fastqc](https://github.com/s-andrews/FastQC)
