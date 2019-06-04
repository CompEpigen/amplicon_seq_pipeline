# amplicon_seq_pipeline
Pipeline to process amplicon sequencing results.

Getting Started
===============

## Depedencies
- [miniconda](https://conda.io/miniconda.html)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [bwa mem](http://bio-bwa.sourceforge.net/)
- [samtools](http://samtools.sourceforge.net/)
- [multiqc](https://multiqc.info/)

## Installation
Clone the pipeline by issuing:
```bash
git clone https://github.com/HeyLifeHD/amplicon_seq_pipeline/
```

## Input
Change the directories of your working directory and reference genome in the config file. Put your fastq.gz files from an amplicon sequencing experiment into a a folder called "01_fastq_raw/" inside your specified working directory. In addition specify the location of the config file in the Snakefile. Now you are ready to go!
