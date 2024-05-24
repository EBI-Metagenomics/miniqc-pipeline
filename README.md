# Microbiome Informatics miniQC pipeline

This pipeline performs quality control (QC) on raw metagenomics reads and decontaminates them using the BWA-MEM2 aligner with a provided reference genome. The pipeline is implemented using Nextflow, enabling scalable and reproducible scientific workflows using software containers.

The code was taken from the [EMBL EBI's Microbiome Informatics team](https://www.ebi.ac.uk/about/teams/microbiome-informatics) [MAGs generation pipeline](https://github.com/EBI-Metagenomics/genomes-generation)

## Introduction

The Metagenomics Reads QC and Decontamination Pipeline is designed to preprocess metagenomic sequencing data by performing quality control and removing contaminant sequences based on a reference genome. The key steps in the pipeline are:
1. Quality control of raw reads using FASTQC.
2. Trimming and filtering of reads using FASTP.
3. Decontamination of reads by aligning to a reference genome using BWA-MEM2.
4. Post-processing quality control using FASTQC.
5. Summarizing the QC results using MultiQC.

## Installation

To run this pipeline, you need to have Nextflow installed. Additionally, Docker or Singularity are required as the pipeline doesn't have Conda support ATM.

### Prerequisites

- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)

## Usage

### Input

The pipeline requires 2 parameters: a samplesheet and a BWA-MEM2 indexed reference genome.

#### Samplesheet

The samplesheet should look like this:

```csv
id,fastq_1,fastq_2
test,/path/to/fastq_1.fastq,
```

It supports paired and single ends. For single ends, leave the third column (fastq_2) empty.

### Execution

```
$ nextflow run ebi-metagenomics/miniqc-pipeline
 N E X T F L O W   ~  version 24.04.1

Launching `main.nf` [hungry_keller] DSL2 - revision: 7a31cd92ec

Typical pipeline command:

  nextflow run ebi-metagenomics/miniqc-pipeline --samplesheet input_file.csv --ref_genome reference.fasta

Input/output options
  --samplesheet                 [string]  Path to comma-separated file containing information about the assemblies, raw reads with the prefix to be 
                                          used. 
  --outdir                      [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud 
                                          infrastructure. [default: results] 

Quality Control
  --merge_pairs                 [boolean] Merge the paired reads after QC

Decontamination
  --ref_genome                  [string]  Reference genome .fasta file, the bwa-mem2 index should be in the same folder

Generic options
  --singularity_cachedir        [string]  null
  --multiqc_title               [string]  Custom title for the MultiQC report.
  --multiqc_methods_description [string]  Custom MultiQC yaml file containing HTML including a methods description.
```