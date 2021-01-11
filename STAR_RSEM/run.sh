#!/bin/bash
mkdir -p /oak/stanford/groups/akundaje/annashch/cd47_RNA/STAR_$1
rm -r /oak/stanford/groups/akundaje/annashch/cd47_RNA/STAR_$1/tmp
mkdir -p /oak/stanford/groups/akundaje/annashch/cd47_RNA/RSEM_$1

#./align.sh /oak/stanford/projects/genomics/data/201027_COOPER_0317_AHCWY3BBXY-SFGF-Chetty-SB-15429/$1 /oak/stanford/groups/akundaje/annashch/cd47_RNA/STAR_$1/$1 /oak/stanford/groups/akundaje/annashch/cd47_RNA/STAR_$1/tmp
./quant.sh /oak/stanford/groups/akundaje/annashch/cd47_RNA/STAR_$1/$1\Aligned.toTranscriptome.out.bam /oak/stanford/groups/akundaje/annashch/cd47_RNA/RSEM_$1/RSEM_$1


