#!/bin/bash

cd /home/alopez.1/XPAND/scRNA/

## Create sample_sheet

## Extract fastq files
cellranger mkfastq --run=260130_VH00461_756_AACYVJVHV --output-dir="Fastq" --samplesheet="sample_RNA.csv"

## Count reads
cellranger count --id="LV12_2dpt" --reference="GRCh38_LV/" --fastqs="Fastq"
