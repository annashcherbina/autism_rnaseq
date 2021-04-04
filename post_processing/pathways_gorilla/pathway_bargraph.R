rm(list=ls())
library(ggplot2)
source("~/helpers.R")
asdn_tdn=read.table("gorilla_asdn_vs_tdn.txt",header=TRUE,sep='\t')
asdn_tdn$FDR=-10*log10(asdn_tdn$FDR)
asdn_tdn$GO=factor(asdn_tdn$GO,levels=asdn_tdn$GO)

asddm_asdn=read.table("gorilla_asddm_vs_tdn.txt",header=TRUE,sep='\t')
asddm_asdn$GO=factor(asddm_asdn$GO,levels=asddm_asdn$GO)
asddm_asdn$FDR=-10*log10(asddm_asdn$FDR)

p1=ggplot(data=asdn_tdn,
          aes(x=GO,
              y=FDR))+
  geom_bar(stat='identity')+ 
    coord_flip()+
    ylab("-10log10(FDR)")+
    xlab("GO Term")+
    ggtitle("Enriched GO Terms\n ASD-N vs TDN")+
    theme_bw(15)
p2=ggplot(data=asddm_asdn,
          aes(x=GO,
              y=FDR))+
  geom_bar(stat='identity')+ 
  coord_flip()+
  ylab("-10log10(FDR)")+
  xlab("GO Term")+
  ggtitle("Enriched GO Terms\n ASD-DM vs TDN")+
  theme_bw(15)
multiplot(p1,p2,cols=2)