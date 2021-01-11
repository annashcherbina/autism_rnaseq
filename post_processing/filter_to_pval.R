rm(list=ls())
library(ggplot2)
source("~/helpers.R")
pthresh=0.000005

npc_asddm_asdn=read.table("diff_genes/NPC_ASDDM_ASDN.diff.tsv",header=TRUE,sep='\t')
npc_asddm_asdn$Gene=row.names(npc_asddm_asdn)
npc_asddm_asdn$sig=npc_asddm_asdn$adj.P.Val<=pthresh & abs(npc_asddm_asdn$logFC)>=1
npc_asddm_asdn$log10Pval=-1*log10(npc_asddm_asdn$adj.P.Val)
npc_asddm_asdn$fcpos=npc_asddm_asdn$logFC>0
up=sum(npc_asddm_asdn$sig==TRUE & npc_asddm_asdn$fcpos==TRUE)
down=sum(npc_asddm_asdn$sig==TRUE & npc_asddm_asdn$fcpos==FALSE)


npc_asddm_tdn=read.table("diff_genes/NPC_ASDDM_TDN.diff.tsv",header=TRUE,sep='\t')
npc_asddm_tdn$Gene=row.names(npc_asddm_tdn)
npc_asddm_tdn$sig=npc_asddm_tdn$adj.P.Val<=pthresh & abs(npc_asddm_tdn$logFC)>=1
npc_asddm_tdn$log10Pval=-1*log10(npc_asddm_tdn$adj.P.Val)
npc_asddm_tdn$fcpos=npc_asddm_tdn$logFC>0
up=sum(npc_asddm_tdn$sig==TRUE & npc_asddm_tdn$fcpos==TRUE)
down=sum(npc_asddm_tdn$sig==TRUE & npc_asddm_tdn$fcpos==FALSE)


npc_asdn_tdn=read.table("diff_genes/NPC_ASDN_TDN.diff.tsv",header=TRUE,sep='\t')
npc_asdn_tdn$Gene=row.names(npc_asdn_tdn)
npc_asdn_tdn$sig=npc_asdn_tdn$adj.P.Val<=pthresh & abs(npc_asdn_tdn$logFC)>=1
npc_asdn_tdn$log10Pval=-1*log10(npc_asdn_tdn$adj.P.Val)
npc_asdn_tdn$fcpos=npc_asdn_tdn$logFC>0
up=sum(npc_asdn_tdn$sig==TRUE & npc_asdn_tdn$fcpos==TRUE)
down=sum(npc_asdn_tdn$sig==TRUE & npc_asdn_tdn$fcpos==FALSE)

#write significant subsets 
write.table(npc_asddm_asdn[npc_asddm_asdn$sig==TRUE,],file=paste("diff_genes/pval.lt.",pthresh,".lfc.gt.1.NPC_ASDDM_ASDN.diff.tsv",sep=""),sep='\t')
write.table(npc_asddm_tdn[npc_asddm_tdn$sig==TRUE,],file=paste("diff_genes/pval.lt.",pthresh,".lfc.gt.1.NPC_ASDDM_TDN.diff.tsv",sep=""),sep='\t')
write.table(npc_asdn_tdn[npc_asdn_tdn$sig==TRUE,],file=paste("diff_genes/pval.lt.",pthresh,".lfc.gt.1.NPC_ASDN_TDN.diff.tsv",sep=""),sep='\t')
