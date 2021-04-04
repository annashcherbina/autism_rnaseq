library(VennDiagram)
asddm_vs_tdn=read.table("NPC_ASDDM_TDN.diff.genes.txt")$V1
asddm_vs_asdn=read.table("NPC_ASDDM_ASDN.diff.genes.txt")$V1
asdn_vs_tdn=read.table("NPC_ASDN_TDN.diff.genes.txt")$V1
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(asddm_vs_tdn, asdn_vs_tdn, asddm_vs_asdn),
  category.names = c("ASDDM_vs_TDN" , "ASDN_vs_TDN" , "ASDDM_vs_ASDN"),
  filename = 'Venn.png',
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