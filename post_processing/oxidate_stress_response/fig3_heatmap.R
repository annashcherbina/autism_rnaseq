library(ggplot2)
library(gplots)
library("gplots")
library("devtools")
library(pheatmap)
data=read.table("oxid.NPC_only.corrected_tpm.txt",header=T,sep='\t',row.names = 1,check.names = F)
batches=read.table("batches.txt",header=TRUE,sep='\t',check.names = F)
batches=batches[batches$Cell=="NPC",]
batches=batches[order(batches$Condition),]

data=data[,batches$TechRep]
annot=as.data.frame(batches$Condition)
colnames(annot)="Condition"
rownames(annot)=batches$TechRep
pheatmap(as.matrix(data), 
         display_numbers = F,
         cluster_rows = T,
         cluster_cols = T,
         scale='row',
         main='Oxidate Stress Gene Expression in NPC samples \n (Row Z-score of Corrected TPM)',
         annotation_col = annot)

# diff genes only 
data=read.table("oxid.de.summary.txt",header=T,sep='\t',row.names = 1,check.names = F)
pheatmap(as.matrix(data), 
         display_numbers = F,
         cluster_rows = F,
         cluster_cols = F,
         scale='none',
         main='Oxidate Stress differential genes, log2FC')
