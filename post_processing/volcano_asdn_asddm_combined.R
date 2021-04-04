rm(list=ls())
library(ggplot2)
source("~/helpers.R")
npc=read.table("diff_genes/NPC.diff.tsv",header=TRUE,sep='\t')
npc$Gene=row.names(npc)
npc$sig=npc$adj.P.Val<=0.01 & abs(npc$logFC)>=1
npc$log10Pval=-1*log10(npc$adj.P.Val)
npc$fcpos=npc$logFC>0
up=sum(npc$sig==TRUE & npc$fcpos==TRUE)
down=sum(npc$sig==TRUE & npc$fcpos==FALSE)
p1=ggplot(npc,
       aes(x=logFC,
           y=log10Pval,
           color=sig))+
  geom_point(size=0.5,alpha=0.5)+
  xlab("logFC")+
  ylab("-log10(adj.p.Val)")+
  scale_color_manual(values=c("#000000","#FF0000"))+
  theme_bw()+
  xlim(-100,100)+
  ggtitle(paste("NPC\nASD vs TDN \nsig=",sum(npc$sig),",\nhigher in ASD=",up,",\nhigher in TDN=",down,sep=''))


ipsc=read.table("diff_genes/IPSC.diff.tsv",header=TRUE,sep='\t')
ipsc$Gene=row.names(ipsc)
ipsc$sig=ipsc$adj.P.Val<=0.01 & abs(ipsc$logFC)>=1
ipsc$log10Pval=-1*log10(ipsc$adj.P.Val)
ipsc$fcpos=ipsc$logFC>0
up=sum(ipsc$sig==TRUE & ipsc$fcpos==TRUE)
down=sum(ipsc$sig==TRUE & ipsc$fcpos==FALSE)
p2=ggplot(ipsc,
          aes(x=logFC,
              y=log10Pval,
              color=sig))+
  geom_point(size=0.5,alpha=0.5)+
  xlab("logFC")+
  ylab("-log10(adj.p.Val)")+
  scale_color_manual(values=c("#000000","#FF0000"))+
  theme_bw()+
  xlim(-100,100)+
  ggtitle(paste("IPSC\nASD vs TDN \nsig=",sum(ipsc$sig),",\nhigher in ASD=",up,",\nhigher in TDN=",down,sep=''))





png("number_of_differential_genes_across_ASD_vs_TDN_comparisons.png",width=12,height=4,units='in',res=300)
multiplot(p1,p2,cols=2)
dev.off()
