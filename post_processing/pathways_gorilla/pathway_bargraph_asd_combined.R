rm(list=ls())
library(ggplot2)
source("~/helpers.R")
npc=read.table("gorilla_npc.txt",header=TRUE,sep='\t')
npc$FDR=-10*log10(npc$FDR)
npc$GO=factor(npc$GO,levels=npc$GO)

p1=ggplot(data=npc,
          aes(x=GO,
              y=FDR))+
  geom_bar(stat='identity')+ 
    coord_flip()+
    ylab("-10log10(FDR)")+
    xlab("GO Term")+
    ggtitle("Enriched GO Terms\n ASD combined vs TDN")+
    theme_bw(15)
