rm(list=ls())
library(ggplot2)
source("~/helpers.R")
pthresh=0.001
lfc_thresh=2

npc_asddm_asdn=read.table("diff_genes/ASDDM_ASDN.diff.tsv",header=TRUE,sep='\t')
npc_asddm_asdn$Gene=row.names(npc_asddm_asdn)
npc_asddm_asdn$sig=npc_asddm_asdn$adj.P.Val<=pthresh & abs(npc_asddm_asdn$logFC)>=lfc_thresh
npc_asddm_asdn$log10Pval=-1*log10(npc_asddm_asdn$adj.P.Val)
npc_asddm_asdn$fcpos=npc_asddm_asdn$logFC>0
up=sum(npc_asddm_asdn$sig==TRUE & npc_asddm_asdn$fcpos==TRUE)
down=sum(npc_asddm_asdn$sig==TRUE & npc_asddm_asdn$fcpos==FALSE)
p1=ggplot(npc_asddm_asdn,
       aes(x=logFC,
           y=log10Pval,
           color=sig))+
  geom_point(size=0.5,alpha=0.5)+
  xlab("logFC")+
  ylab("-log10(adj.p.Val)")+
  scale_color_manual(values=c("#000000","#FF0000"))+
  theme_bw()+
  xlim(-100,100)+
  ggtitle(paste("NPC\nASD-DM vs ASD-N\nsig=",sum(npc_asddm_asdn$sig),",\nhigher in ASD-DM=",up,",\nhigher in ASD-N=",down,sep=''))


npc_asddm_tdn=read.table("diff_genes/ASDDM_TDN.diff.tsv",header=TRUE,sep='\t')
npc_asddm_tdn$Gene=row.names(npc_asddm_tdn)
npc_asddm_tdn$sig=npc_asddm_tdn$adj.P.Val<=pthresh & abs(npc_asddm_tdn$logFC)>=lfc_thresh
npc_asddm_tdn$log10Pval=-1*log10(npc_asddm_tdn$adj.P.Val)
npc_asddm_tdn$fcpos=npc_asddm_tdn$logFC>0
up=sum(npc_asddm_tdn$sig==TRUE & npc_asddm_tdn$fcpos==TRUE)
down=sum(npc_asddm_tdn$sig==TRUE & npc_asddm_tdn$fcpos==FALSE)
p2=ggplot(npc_asddm_tdn,
          aes(x=logFC,
              y=log10Pval,
              color=sig))+
  geom_point(size=0.5,alpha=0.5)+
  xlab("logFC")+
  ylab("-log10(adj.p.Val)")+
  scale_color_manual(values=c("#000000","#FF0000"))+
  theme_bw()+
  xlim(-100,100)+
  ggtitle(paste("NPC\nASD-DM vs TDN\nsig=",sum(npc_asddm_tdn$sig),",\nhigher in ASD-DM=",up,",\nhigher in TDN=",down,sep=''))


npc_asdn_tdn=read.table("diff_genes/ASDN_TDN.diff.tsv",header=TRUE,sep='\t')
npc_asdn_tdn$Gene=row.names(npc_asdn_tdn)
npc_asdn_tdn$sig=npc_asdn_tdn$adj.P.Val<=pthresh & abs(npc_asdn_tdn$logFC)>=lfc_thresh
npc_asdn_tdn$log10Pval=-1*log10(npc_asdn_tdn$adj.P.Val)
npc_asdn_tdn$fcpos=npc_asdn_tdn$logFC>0
up=sum(npc_asdn_tdn$sig==TRUE & npc_asdn_tdn$fcpos==TRUE)
down=sum(npc_asdn_tdn$sig==TRUE & npc_asdn_tdn$fcpos==FALSE)
p3=ggplot(npc_asdn_tdn,
          aes(x=logFC,
              y=log10Pval,
              color=sig))+
  geom_point(size=0.5,alpha=0.5)+
  xlab("logFC")+
  ylab("-log10(adj.p.Val)")+
  scale_color_manual(values=c("#000000","#FF0000"))+
  theme_bw()+
  xlim(-100,100)+
  ggtitle(paste("NPC\nASD-N vs TDN\nsig=",sum(npc_asdn_tdn$sig),",\nhigher in ASD-N=",up,",\nhigher in TDN=",down,sep=''))

png("number_of_differential_genes_across_NPC_comparisons.png",width=12,height=4,units='in',res=300)
multiplot(p3,p2,p1,cols=3)
dev.off()

#write significant subsets 
write.table(npc_asddm_asdn[npc_asddm_asdn$sig==TRUE,],file=paste("diff_genes/pval.lt.",pthresh,".lfc.gt.",lfc_thresh,".NPC_ASDDM_ASDN.diff.tsv",sep=""),sep='\t')
write.table(npc_asddm_tdn[npc_asddm_tdn$sig==TRUE,],file=paste("diff_genes/pval.lt.",pthresh,".lfc.gt.",lfc_thresh,".NPC_ASDDM_TDN.diff.tsv",sep=""),sep='\t')
write.table(npc_asdn_tdn[npc_asdn_tdn$sig==TRUE,],file=paste("diff_genes/pval.lt.",pthresh,".lfc.gt.",lfc_thresh,".NPC_ASDN_TDN.diff.tsv",sep=""),sep='\t')
