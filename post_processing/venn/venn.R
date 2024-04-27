library(VennDiagram)
#asddm_vs_tdn=row.names(read.table("../pval.lt.0.001.lfc.gt.2.NPC_ASDDM_TDN.diff.tsv"))
#asddm_vs_asdn=row.names(read.table("../pval.lt.0.001.lfc.gt.2.NPC_ASDDM_ASDN.diff.tsv"))
#asdn_vs_tdn=row.names(read.table("../pval.lt.0.001.lfc.gt.2.NPC_ASDN_TDN.diff.tsv"))
asddm_vs_tdn=row.names(read.table("../diff_genes/ConditionASDDM_vs_ConditionTDN.de.sig.tsv"))
asddm_vs_asdn=row.names(read.table("../diff_genes/ConditionASDDM_vs_ConditionASDN.de.sig.tsv"))
asdn_vs_tdn=row.names(read.table("../diff_genes/ConditionASDN_vs_ConditionTDN.de.sig.tsv"))

library(RColorBrewer)
myCol <- c('#377eb8','#4daf4a','#e41a1c')
venn.diagram(
  x = list(asddm_vs_tdn, asdn_vs_tdn, asddm_vs_asdn),
  category.names = c("ASDDM_vs_TDN" , "ASDN_vs_TDN" , "ASDDM_vs_ASDN"),
  filename = 'Venn.svg',
  output=TRUE,
  # Output features
  imagetype="svg" ,
  height = 800 , 
  width = 800 , 
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)