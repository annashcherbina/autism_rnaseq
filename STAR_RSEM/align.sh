#!/bin/bash
/home/groups/akundaje/annashch/bulk-rna-seq/rna_seq_pipeline/map_rna_star.sh --R1 $1\_R1_001.fastq.gz --R2 $1\_R2_001.fastq.gz --starGenomeDir /oak/stanford/groups/akundaje/refs/STAR/GRCh38 --ncpus 16 --outputbase $2 --outTmpDir $3


