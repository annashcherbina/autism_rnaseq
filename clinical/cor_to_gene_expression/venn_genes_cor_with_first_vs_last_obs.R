library(VennDiagram)
asddm_vs_tdn_first=row.names(read.table('spearman.ConditionASDDM_vs_ConditionTDN.first.tsv',header=T,sep='\t'))
asddm_vs_tdn_last=row.names(read.table('spearman.ConditionASDDM_vs_ConditionTDN.last.tsv',header=T,sep='\t'))

asdn_vs_tdn_first=row.names(read.table('spearman.ConditionASDN_vs_ConditionTDN.first.tsv',header=T,sep='\t'))
asdn_vs_tdn_last=row.names(read.table('spearman.ConditionASDN_vs_ConditionTDN.last.tsv',header=T,sep='\t'))


asddm_vs_asdn_first=row.names(read.table('spearman.ConditionASDDM_vs_ConditionASDN.first.tsv',header=T,sep='\t'))
asddm_vs_asdn_last=row.names(read.table('spearman.ConditionASDDM_vs_ConditionASDN.last.tsv',header=T,sep='\t'))

library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel2")
venn.diagram(
  x = list(asddm_vs_tdn_first, asddm_vs_tdn_last),
  category.names = c("ASDDM_vs_TDN_first","ASDDM_vs_TDN_last"),
  filename = 'Venn.ASDDM_vs_TDN_first_vs_last.png',
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 800 ,
  width = 800 ,
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055)
)