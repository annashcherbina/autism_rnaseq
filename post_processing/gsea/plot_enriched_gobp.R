library(pheatmap)
data=read.table("filtered/merged.filtered.0.01.0.5.fdr_thresholded_nes.tsv",header=T,sep='\t',row.names = 1)
colnames(data)=c("ASDDM_vs_ASDN","ASDDM_vs_TDN","ASDN_vs_TDN")
data[is.na(data)]=0
pheatmap(t(data), 
         display_numbers = T,
         cluster_rows = T,
         cluster_cols = T,
         main='Normalized enrichmed score for significant GO:BP (p.adj < 0.01) ')