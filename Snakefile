#Snakmake Pipeline


import os
from os import path

configfile: "/home/epicwl/c010-datasets/Internal/GCTB/MiSeq/raw_data/DKFZ_run/config.yml"
workdir: config["workdir_top"]

import sys

def message(mes):
  sys.stderr.write("|--- " + mes + "\n")


import glob
import re

SAMPLES, = glob_wildcards("01_fastq_raw/{sam}_R1.fastq.gz")

NB_SAMPLES = len(SAMPLES)

message(str(NB_SAMPLES))

for sample in SAMPLES:
  message("Sample " + str(sample) + " will be processed")

rule all:
    input:
        expand("03_bwa/{sample}_sorted.sam", sample=SAMPLES)

rule multiQC:
    message:
        "Perform MultiQC"
    shell:"""
        multiqc --title "Raw FastQC reports" --filename Raw_FastQC_report --outdir MultiQC/ ./
        """

rule trim:
    message:
        "Preprocessing of fastq files"
    input:
        fwd = "01_fastq_raw/{sample}_R1.fastq.gz",
        rev = "01_fastq_raw/{sample}_R2.fastq.gz"
    output:
        fwd = "02_fastq_trimmed/{sample}_R1_val_1.fq.gz",
        rev = "02_fastq_trimmed/{sample}_R2_val_2.fq.gz"
    shell:"""
        trim_galore --paired --nextera  {input.fwd} {input.rev} -o 02_fastq_trimmed/
        touch {output.fwd}
        touch {output.rev}
        """


rule align:
#maybe path to bowtie still has to be added
    message:
        "Align trimmed fastq files"
    input:
        fwd = "02_fastq_trimmed/{sample}_R1_val_1.fq.gz",
        rev = "02_fastq_trimmed/{sample}_R2_val_2.fq.gz"
    output:
        sam =  "03_bwa/{sample}.sam"
    params:
        ref= config["BWA_reference"]
    shell:"""
        bwa mem {params.ref} -t {threads} {input.fwd} {input.rev} > {output.sam}
        #samtools sort -@ 4 -o {output.sam} {output.sam}
        """

rule sort:
    message:
        "Sorting sam files"
    input:
        sam =  "03_bwa/{sample}.sam"
    output:
        sam =  "03_bwa/{sample}_sorted.sam"
    shell:"""
        samtools sort -@ 4 -o {output.sam} {input.sam}
        """