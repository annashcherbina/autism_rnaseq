#sort by nes from most positive to most negative 
#go 
go_asddm_asdn=read.table("plot_subset/ASDDM_ASDN.diff.tsv.GSEA.tsv.sig.GOBP.tsv",header=T,sep='\t')
go_asddm_asdn$Term=factor(go_asddm_asdn$Term,levels=go_asddm_asdn$Term)

go_asddm_tdn=read.table("plot_subset/ASDDM_TDN.diff.tsv.GSEA.tsv.sig.GOBP.tsv",header=T,sep='\t')
go_asddm_tdn$Term=factor(go_asddm_tdn$Term,levels=go_asddm_tdn$Term)

go_asdn_tdn=read.table('plot_subset/ASDN_TDN.diff.tsv.GSEA.tsv.sig.GOBP.tsv',header=T,sep='\t')
go_asdn_tdn$Term=factor(go_asdn_tdn$Term,levels=go_asdn_tdn$Term)


library(ggplot2)

p1=ggplot(go_asddm_asdn,
          aes(x=Term,
              y=nes))+
  geom_bar(stat='identity')+
  xlab('Term')+
  ylab("normalized enrichment score")+
  coord_flip()+
  theme(text = element_text(size=rel(3)),
        plot.title = element_text(size = 10, face = "bold"))+
  ggtitle("GO BP, ASDDM vs ASDN")

p2=ggplot(go_asddm_tdn,
          aes(x=Term,
              y=nes))+
  geom_bar(stat='identity')+
  xlab('Term')+
  ylab("normalized enrichment score")+
  coord_flip()+
  theme(text = element_text(size=rel(3)),
        plot.title = element_text(size = 10, face = "bold"))+
  ggtitle("GO BP, ASDDM vs TDN")

p3=ggplot(go_asdn_tdn,
          aes(x=Term,
              y=nes))+
  geom_bar(stat='identity')+
  xlab('Term')+
  ylab("normalized enrichment score")+
  coord_flip()+
  theme(text = element_text(size=rel(3)),
        plot.title = element_text(size = 10, face = "bold"))+
  ggtitle("GO BP, ASDN vs TDN")

multiplot(p1,p2,p3,cols=3)