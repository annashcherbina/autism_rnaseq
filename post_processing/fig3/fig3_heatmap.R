library(ggplot2)
library(gplots)
library("gplots")
library("devtools")
fig3_genes=read.table("fig3_genes.txt",header=F,sep='\t')$V1
data=read.table("../NPC_only.corrected_tpm.txt",header=T,sep='\t',row.names = 1)
common=intersect(rownames(data),fig3_genes)
data=data[common,]
batches=read.table("../merged_rsem/batches.txt",header=TRUE,sep='\t',check.names = F)
batches=batches[batches$Cell=="NPC",]
batches=batches[order(batches$Condition),]

data=data[,batches$TechRep]
pheatmap(as.matrix(data), 
         display_numbers = F,
         cluster_rows = T,
         cluster_cols = F,
         scale='row',
         main='row Z in NPC for genes from Palmer et al.',
         annotation_col = as.data.frame(batches$Condition))

