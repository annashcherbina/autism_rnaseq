library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel2")
for(comparison in c("ASDN_TDN","ASDDM_TDN","ASDDM_ASDN"))
{
set1=row.names(read.table(paste0('spearman.',comparison,'.first.tsv'),header=T,sep='\t'))
set2=row.names(read.table(paste0('spearman.',comparison,'.last.tsv'),header=T,sep='\t'))
venn.diagram(
  x = list(set1, set2),
  category.names = c(paste(comparison,'first'),
                     paste(comparison,'last')),
  filename = paste0(comparison,'_venn_diagramm.png'),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2]
)
}