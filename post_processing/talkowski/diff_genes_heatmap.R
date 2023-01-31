library(ggplot2)
library(gplots)
library("gplots")
library("devtools")

asdn_vs_tdn=read.table('ConditionASDN_vs_ConditionTDN.de.sig.tsv',header=T,sep='\t',check.names = F)
asdn_vs_tdn=asdn_vs_tdn[order(asdn_vs_tdn$adj.P.Val),]
asddm_vs_tdn=read.table('ConditionASDDM_vs_ConditionTDN.de.sig.tsv',header=T,sep='\t',check.names = F)
asddm_vs_tdn=asddm_vs_tdn[order(asddm_vs_tdn$adj.P.Val),]

asddm_vs_asdn=read.table('ConditionASDDM_vs_ConditionASDN.de.sig.tsv',header=T,sep='\t',check.names=F)
asddm_vs_asdn=asddm_vs_asdn[order(asddm_vs_asdn$adj.P.Val),]

top_sig_genes=unique(c(asddm_vs_tdn$Gene,asddm_vs_asdn$Gene, asdn_vs_tdn$Gene))
data=read.table("../NPC_only.corrected_tpm.txt",header=T,sep='\t',row.names = 1, check.names = F)[top_sig_genes,]

batches=read.table("../merged_rsem/batches.txt",header=TRUE,sep='\t',check.names = F)
batches=batches[batches$Cell=="NPC",]
batches=batches[order(batches$Condition),]
library(pheatmap)
data=data[,batches$TechRep]

#filter to talkowski
talkowski_genes=read.table("talkowski_genes.txt",header=F)$V1
common=intersect(talkowski_genes,rownames(data))
data=data[common,]
annotation_row=subset(batches,select=c("Condition","TechRep"))
rownames(annotation_row)=annotation_row$TechRep
annotation_row$TechRep=NULL

pheatmap(as.matrix(t(data)), 
         display_numbers = F,
         cluster_rows = T,
         cluster_cols = T,
         main='corrected TPM DE genes shared w/ Talkowski et al',
         annotation_row = annotation_row)

