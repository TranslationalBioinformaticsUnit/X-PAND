#!/bin/sh
## Concatenate virus fasta files with human genome fasta
cat Xpand_Ref/LV.fa Xpand_Ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa > Xpand_Ref/Homo_sapiens.GRCh38.dna.primary_assembly_LV.fa
cat Xpand_Ref/RV.fa Xpand_Ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa > Xpand_Ref/Homo_sapiens.GRCh38.dna.primary_assembly_RV.fa

## Add gene information to the GTF files

## Create config file

cellranger-arc mkref --config=Xpand_Ref/GRCh38_LV.config
cellranger-arc mkref --config=Xpand_Ref/GRCh38_RV.config
