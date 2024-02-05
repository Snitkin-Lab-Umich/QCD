# QCD - Quality Control and Contamination Detection workflow.

QCD is a smakemake worflow for microbial illumina sequencing quality control and contamination detection. As a SOP part of Snitkin lab, this pipeline should be run on raw sequencing data as soon as the data is available from the sequencing core department. Apart from QCing raw sequencing data, it performs numerous downstream tasks such as AMR gene detection, SPAdes genome assembly, MLST detection and Assembly annotation. 

## Installation

> Clone the github directory onto your system.

```
git clone https://github.com/Snitkin-Lab-Umich/QCD.git

```
> Load snakemake module from Great Lakes modules

```
module load snakemake
```

> Customise snakemake configuration settings in config/config.yaml file as per your needs and create a sample list file for your project - sample.tsv


## Quick start

### Run QCD on a set of samples.

Run QCD locally

```
snakemake -s QCD.smk -p --configfile config/config.yaml --cores all
```

Run QCD on Great lakes HPC

```
snakemake -s QCD.smk -p --use-conda -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes}  -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 10000
```

![Alt text](./QCD_dag.svg)

### Gather Summary files and generate a report. 
```
snakemake -s QCD_report.smk -p --configfile config/config.yaml --cores all
```
![Alt text](./QCD_report_dag.svg)

## Dependencies

### Near Essential
* [Snakemake>=7.32.4](https://snakemake.readthedocs.io/en/stable/#)
* [Conda](https://docs.conda.io/en/latest/)

All the necessary software stack required for the workflow will be installed using conda package manager.

### Tool stack used in workflow

* [fastq-scan](https://github.com/rpetit3/fastq-scan)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [SPades](https://github.com/ablab/spades)
* [AMRFinderPlus](https://github.com/ncbi/amr)
* [bioawk](https://github.com/lh3/bioawk)
* [Prokka](https://github.com/tseemann/prokka)
* [mlst](https://github.com/tseemann/mlst)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
