#visual inspection to remove differential genes with var / mean (max across conditions) > 0.5 
rm(list=ls())
library(ggplot2)
var_thresh=0.75
data=read.table("NPC_only.corrected_tpm.txt",header=TRUE,sep='\t',row.names = 1,check.names=FALSE)
batches=read.table("merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)
conditions = c("TDN","ASDN","ASDDM")
for(comparison in c("ConditionASDDM_vs_ConditionASDN","ConditionASDDM_vs_ConditionTDN","ConditionASDN_vs_ConditionTDN"))
{
  de=read.table(paste0("diff_genes/",comparison,".de.sig.tsv"),header=T,sep='\t')
  de$max_var_mean_ratio=pmax(de$TDN_var_mean_ratio,
                             de$ASDDM_var_mean_ratio,
                             de$ASDN_var_mean_ratio)
  genes_of_interest=rownames(de)[de$max_var_mean_ratio>var_thresh]
  #subset data to genes of interest
  sub_data=data[genes_of_interest,]
  sub_data=as.data.frame(t(sub_data))
  sub_data=merge(sub_data,batches,by=0)
  sub_data$Condition=factor(sub_data$Condition, levels=conditions)
  plot.list = lapply(genes_of_interest, function(gene){
  p <- ggplot(data=sub_data,
              aes(x=sub_data$Condition,
                  y=sub_data[[gene]]))+
    geom_boxplot(color='#FF0000',outlier.shape=NA)+
    geom_jitter()+
    ylab("TPM")+
    xlab("Condition")+
    ggtitle(gene)+
    theme_bw(5)
  return(p)
})
library(gridExtra)
grid.arrange(grobs = plot.list , ncol=10,nrow=6)
}