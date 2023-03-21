rm(list=ls())
library(stringr)
library(ggplot2)

wrap_size=50
#ASDN vs TDN
asdn_tdn=read.table("ConditionASDN_vs_ConditionTDN.de.tsv.GSEA.tsv",header=T,sep='\t')
asdn_tdn$Term=str_wrap(str_replace(str_replace(tolower(asdn_tdn$Term),'gobp_','GOBP:'),'_',' '),wrap_size)
asdn_tdn$log_fdr=-1*log10(asdn_tdn$fdr+1e-10)
asdn_tdn$Comparison="ASD-N vs TD-N"
asdn_tdn$Term=as.character(asdn_tdn$Term)


#ASDDM vs ASDN
asddm_asdn=read.table("ConditionASDDM_vs_ConditionASDN.de.tsv.GSEA.tsv",header=T,sep='\t')
asddm_asdn$Term=str_wrap(str_replace(str_replace(tolower(asddm_asdn$Term),'gobp_','GOBP:'),'_',' '),wrap_size)
asddm_asdn$log_fdr=-1*log10(asddm_asdn$fdr+1e-10)
asddm_asdn$Comparison="ASD-DM vs ASD-N"
asddm_asdn$Term=as.character(asddm_asdn$Term)

#ASDDM vs TDN
asddm_tdn=read.table("ConditionASDDM_vs_ConditionTDN.de.tsv.GSEA.tsv",header=T,sep='\t')
asddm_tdn$Term=str_wrap(str_replace(str_replace(tolower(asddm_tdn$Term),'gobp_','GOBP:'),'_',' '),wrap_size)
asddm_tdn$log_fdr=-1*log10(asddm_tdn$fdr+1e-10)
asddm_tdn$Term=as.character(asddm_tdn$Term)
asddm_tdn$Comparison="ASD-DM vs TD-N"
asddm_tdn$Term=as.character(asddm_tdn$Term)


merged=rbind(asddm_tdn,asdn_tdn,asddm_asdn)
merged=merged[order(merged$nes,decreasing = T),]
library(stringr)
merged$Term=str_replace_all(merged$Term,'GOBP:','')
merged$Term=str_replace_all(merged$Term,'_',' ')
merged$Term=str_to_title(merged$Term)


up=merged[merged$nes > 0,]
down=merged[merged$nes<0,]

up$Term=factor(up$Term,levels=unique(up$Term))

p1=ggplot(up,
       aes(x=Term,
           y=nes,
           group=Comparison,
           fill=Comparison))+
  geom_bar(stat='identity',position='dodge')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("Enriched GO Term")+ylab("Normalized Effect Size (GSEA)")+
  scale_fill_manual(name="-log10(FDR)",values=c('#377eb8','#4daf4a','#e41a1c'))+
  ggtitle("GSEA GO Term Enrichments:\nUpregulated in ASD relative to TD") + 
  theme_bw(16)+
  theme(axis.text = element_text(face="bold"))+
  theme(legend.position="bottom")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

down$Term=factor(down$Term,levels=unique(down$Term))

p2=ggplot(down,
       aes(x=Term,
           y=nes,
           grodown=Comparison,
           fill=Comparison))+
  geom_bar(stat='identity',position='dodge')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("Enriched GO Term")+ylab("Normalized Effect Size (GSEA)")+
  scale_fill_manual(name="-log10(FDR)",values=c('#377eb8','#4daf4a','#e41a1c'))+
  ggtitle("GSEA GO Term Enrichments:\nDownregulated in ASD Relative to TD") + 
  theme_bw(16)+
  theme(axis.text = element_text(face="bold"))+
  theme(legend.position="bottom")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 70))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
library(scater)
multiplot(p1,p2,cols=1)


