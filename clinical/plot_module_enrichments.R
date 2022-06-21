library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
data=read.table('module.enrichments.clinical.assoc.genes.tsv',header=T,sep='\t')
data=data[data$P.val<0.05,]
comparisons=c("ASDN_TDN","ASDDM_TDN","ASDDM_ASDN")
for(comparison in comparisons)
{
  data_comparison=data[data$Comparison==comparison,c('ClinFeat','Module','P.val')]
  data_comparison=dcast(data_comparison,ClinFeat~Module)
  rownames(data_comparison)=data_comparison$ClinFeat
  data_comparison$ClinFeat=NULL
  data_comparison=-1*log10(data_comparison)
  print(Heatmap(as.matrix(data_comparison),
                cluster_rows = F, 
                cluster_columns = F,
                name=paste(comparison,'-log10(p)',sep='\n'),
                row_title = 'Clinical Feature',
                column_title = 'Enriched WGCNA Module'))
}
