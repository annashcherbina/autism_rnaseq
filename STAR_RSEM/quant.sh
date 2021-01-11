#!/bin/bash
echo $1 
echo $2
/home/groups/akundaje/annashch/bulk-rna-seq/rna_seq_pipeline/quant_rna_rsem.sh --rsemGenomeDir /oak/stanford/groups/akundaje/refs/RSEM/GRCh38 --transcriptBAM $1 --ncpus 16 --output $2
