#!/bin/bash

## Create sample_sheet

cellranger-atac mkfastq --run=260205_VH00461_759_222HC3HNX --output-dir="Fastq" --samplesheet="sample_ATAC.csv"

cellranger-atac count --id="RV12_2dpt" --sample="RV12_2dpt" --reference="GRCh38_RV" --fastqs="Fastq" 
