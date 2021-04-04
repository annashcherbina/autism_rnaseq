rm(list=ls())
library(ggplot2)
source("~/helpers.R")
pthresh=0.01

npc=read.table("diff_genes/NPC.diff.tsv",header=TRUE,sep='\t')
npc$Gene=row.names(npc)
npc$sig=npc$adj.P.Val<=pthresh & abs(npc$logFC)>=1
npc$log10Pval=-1*log10(npc$adj.P.Val)
npc$fcpos=npc$logFC>0
up=sum(npc$sig==TRUE & npc$fcpos==TRUE)
down=sum(npc$sig==TRUE & npc$fcpos==FALSE)

#write significant subsets 
write.table(npc[npc$sig==TRUE,],file=paste("diff_genes/pval.lt.",pthresh,".lfc.gt.1.NPC.diff.tsv",sep=""),sep='\t')


ipsc=read.table("diff_genes/IPSC.diff.tsv",header=TRUE,sep='\t')
ipsc$Gene=row.names(ipsc)
ipsc$sig=ipsc$adj.P.Val<=pthresh & abs(ipsc$logFC)>=1
ipsc$log10Pval=-1*log10(ipsc$adj.P.Val)
ipsc$fcpos=ipsc$logFC>0
up=sum(ipsc$sig==TRUE & ipsc$fcpos==TRUE)
down=sum(ipsc$sig==TRUE & ipsc$fcpos==FALSE)

#write significant subsets 
write.table(ipsc[ipsc$sig==TRUE,],file=paste("diff_genes/pval.lt.",pthresh,".lfc.gt.1.ipsc.diff.tsv",sep=""),sep='\t')
