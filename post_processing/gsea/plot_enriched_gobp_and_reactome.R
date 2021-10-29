#sort by nes from most positive to most negative 
#go 
go_asddm_asdn=read.table("ASDDM_ASDN.diff.tsv.GSEA.tsv.sig.GOBP.tsv",header=T,sep='\t')
go_asddm_asdn=go_asddm_asdn[order(go_asddm_asdn$nes,decreasing=T),]
go_asddm_asdn$Term=factor(go_asddm_asdn$Term,levels=go_asddm_asdn$Term)

go_asddm_tdn=read.table("ASDDM_TDN.diff.tsv.GSEA.tsv.sig.GOBP.tsv",header=T,sep='\t')
go_asddm_tdn=go_asddm_tdn[order(go_asddm_tdn$nes,decreasing=T),]
go_asddm_tdn$Term=factor(go_asddm_tdn$Term,levels=go_asddm_tdn$Term)

go_asdn_tdn=read.table('ASDN_TDN.diff.tsv.GSEA.tsv.sig.GOBP.tsv',header=T,sep='\t')
go_asdn_tdn=go_asdn_tdn[order(go_asdn_tdn$nes,decreasing=T),]
go_asdn_tdn$Term=factor(go_asdn_tdn$Term,levels=go_asdn_tdn$Term)

#reactome 
reactome_asddm_asdn=read.table("ASDDM_ASDN.diff.tsv.GSEA.tsv.sig.REACTOME.tsv",header=T,sep='\t')
reactome_asddm_asdn=reactome_asddm_asdn[order(reactome_asddm_asdn$nes,decreasing=T),]
reactome_asddm_asdn$Term=factor(reactome_asddm_asdn$Term,levels=reactome_asddm_asdn$Term)

reactome_asddm_tdn=read.table("ASDDM_TDN.diff.tsv.GSEA.tsv.sig.REACTOME.tsv",header=T,sep='\t')
reactome_asddm_tdn=reactome_asddm_tdn[order(reactome_asddm_tdn$nes,decreasing=T),]
reactome_asddm_tdn$Term=factor(reactome_asddm_tdn$Term,levels=reactome_asddm_tdn$Term)

reactome_asdn_tdn=read.table('ASDN_TDN.diff.tsv.GSEA.tsv.sig.REACTOME.tsv',header=T,sep='\t')
reactome_asdn_tdn=reactome_asdn_tdn[order(reactome_asdn_tdn$nes,decreasing=T),]
reactome_asdn_tdn$Term=factor(reactome_asdn_tdn$Term,levels=reactome_asdn_tdn$Term)


library(ggplot2)

p1=ggplot(go_asddm_asdn,
          aes(x=Term,
              y=nes))+
  geom_bar(stat='identity')+
  xlab('Term')+
  ylab("normalized enrichment score")+
  coord_flip()+
  theme(text = element_text(size=rel(2.5)))+
  ggtitle("GO BP, ASDDM vs ASDN")

p2=ggplot(go_asddm_tdn,
          aes(x=Term,
              y=nes))+
  geom_bar(stat='identity')+
  xlab('Term')+
  ylab("normalized enrichment score")+
  coord_flip()+
  theme(text = element_text(size=rel(2.5)))+
  ggtitle("GO BP, ASDDM vs TDN")

p3=ggplot(go_asdn_tdn,
           aes(x=Term,
               y=nes))+
  geom_bar(stat='identity')+
  xlab('Term')+
  ylab("normalized enrichment score")+
  coord_flip()+
  theme(text = element_text(size=rel(2.5)))+
  ggtitle("GO BP, ASDN vs TDN")


p4=ggplot(reactome_asddm_asdn,
          aes(x=Term,
              y=nes))+
  geom_bar(stat='identity')+
  xlab('Term')+
  ylab("normalized enrichment score")+
  coord_flip()+
  theme(text = element_text(size=rel(2.5)))+
  ggtitle("Reactome BP, ASDDM vs ASDN")

p5=ggplot(reactome_asddm_tdn,
          aes(x=Term,
              y=nes))+
  geom_bar(stat='identity')+
  xlab('Term')+
  ylab("normalized enrichment score")+
  coord_flip()+
  theme(text = element_text(size=rel(2.5)))+
  ggtitle("Reactome BP, ASDDM vs TDN")

p6=ggplot(reactome_asdn_tdn,
          aes(x=Term,
              y=nes))+
  geom_bar(stat='identity')+
  xlab('Term')+
  ylab("normalized enrichment score")+
  coord_flip()+
  theme(text = element_text(size=rel(2.5)))+
  ggtitle("Reactome BP, ASDN vs TDN")


