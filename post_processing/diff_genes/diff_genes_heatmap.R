library(ggplot2)
library(gplots)
library("gplots")
library("devtools")

topn=50
asdn_vs_tdn=read.table('ConditionASDN_vs_ConditionTDN.de.sig.varfilter.tsv',header=T,sep='\t')
asdn_vs_tdn=asdn_vs_tdn[order(asdn_vs_tdn$adj.P.Val),][1:topn,]
asddm_vs_tdn=read.table('ConditionASDDM_vs_ConditionTDN.de.sig.varfilter.tsv',header=T,sep='\t')
asddm_vs_tdn=asddm_vs_tdn[order(asddm_vs_tdn$adj.P.Val),][1:topn,]

asddm_vs_asdn=read.table('ConditionASDDM_vs_ConditionASDN.de.sig.varfilter.tsv',header=T,sep='\t')
asddm_vs_asdn=asddm_vs_asdn[order(asddm_vs_asdn$adj.P.Val),][1:topn,]

top_sig_genes=unique(c(asddm_vs_tdn$Gene,asddm_vs_asdn$Gene, asdn_vs_tdn$Gene))
data=read.table("../NPC_only.corrected_tpm.txt",header=T,sep='\t',row.names = 1)[top_sig_genes,]

batches=read.table("../merged_rsem/batches.txt",header=TRUE,sep='\t',check.names = F)
batches=batches[batches$Cell=="NPC",]
batches=batches[order(batches$Condition),]

data=data[,batches$TechRep]
pheatmap(as.matrix(t(log2(data))), 
         display_numbers = F,
         cluster_rows = F,
         cluster_cols = T,
         main='log2(corrected TPM) in NPC for top 50 sig genes per comparison',
         annotation_row = as.data.frame(batches$Condition))

