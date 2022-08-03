rm(list=ls())
library(EnhancedVolcano)
library(BiocParallel)
source('volcano_wrapper.R')

samples=c("dmso_v_pbs_3","dmso_v_pbs_6","dmso_v_pbs_24","dmso_v_pbs_48",
          "antitgfb_v_pbs_3","antitgfb_v_pbs_6","antitgfb_v_pbs_24","antitgfb_v_pbs_48",
          "tgfb_v_pbs_3","tgfb_v_pbs_6","tgfb_v_pbs_24","tgfb_v_pbs_48",
          "tnfa_v_pbs_3","tnfa_v_pbs_6","tnfa_v_pbs_24","tnfa_v_pbs_48",
          "tgfbtnfa_v_pbs_3","tgfbtnfa_v_pbs_6","tgfbtnfa_v_pbs_24","tgfbtnfa_v_pbs_48")
grobs=list()
for(sample in samples)
  
{
  data=read.table(paste(sample,".diff.tsv",sep=""),header=T,sep='\t',row.names=1)
  data$feature=row.names(data)
  grobs[[sample]]=volcano_wrapper(data,title=sample,pval_thresh=0.05,lfc_thresh=0.1)
}
png(
  "kinetics.25million.volcano.png",
  width = 20,
  height = 30,
  units='in',
  res=300
)
gridExtra::grid.arrange(grobs = grobs, nrow = 5)
dev.off()