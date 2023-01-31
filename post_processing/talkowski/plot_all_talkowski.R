library(ggplot2)
library(gplots)
library("gplots")
library("devtools")
library(tidyr)
data=read.table("../merged_rsem/rna.counts.txt.tpm",header=TRUE,check.names = F)
data=separate(data = data, col = "GeneID", into = c("GeneID", "GeneIDMinor"), sep = "\\.")
talkowski_genes=read.table("talkowski_genes.txt",header=F)$V2
common=intersect(talkowski_genes,data$GeneID)
data=data[data$GeneID %in% common,]
batches=read.table("../merged_rsem/batches.txt",header=TRUE)
rownames(data)=data$GeneName
data$GeneID=NULL
data$GeneIDMinor=NULL
data$GeneName=NULL
batches=batches[order(batches$Cell,batches$Condition),]
data=data[,batches$TechRep]
annotation_row=subset(batches,select=c("Condition","TechRep","Cell"))
rownames(annotation_row)=annotation_row$TechRep
annotation_row$TechRep=NULL

pheatmap(as.matrix(t(log2(data+0.001))), 
         display_numbers = F,
         cluster_rows = F,
         cluster_cols = T,
         main='log2(TPM) all genes from Talkowski et al',
         annotation_row = annotation_row)

#filter to talkowski
annotation_row=subset(batches,select=c("Condition","TechRep"))
rownames(annotation_row)=annotation_row$TechRep
annotation_row$TechRep=NULL

pheatmap(as.matrix(t(data)), 
         display_numbers = F,
         cluster_rows = T,
         cluster_cols = T,
         main='corrected TPM DE genes shared w/ Talkowski et al',
         annotation_row = annotation_row)

