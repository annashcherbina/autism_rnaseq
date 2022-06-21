rm(list=ls())
library(limma)
library(sva)
library(ggplot2)
library(gProfileR)
source('../volcano_wrapper.R')
data=read.table("../post_processing/corrected_tpm.txt",header=TRUE,check.names = F)
batches=read.table("../post_processing/merged_rsem/batches.txt",header=TRUE)
batches$Sample=factor(batches$Sample)

#LIMMA REGRESSION ANALYSIS
first=as.data.frame(scale(read.table('first_vals.txt',header=T,sep='\t',row.names=1)))
first$ACE_MIND_ID=row.names(first)
last=as.data.frame(scale(read.table('last_vals.txt',header=T,sep='\t',row.names=1)))
last$ACE_MIND_ID=row.names(last)
mapfile=read.table('APP iPSC pilot study_ID numbers.tsv',header=T,sep='\t')
batches_anno=merge(batches,mapfile,by.y ='iPSC_ID',by.x='Sample')
batches_first=merge(first,batches_anno,by='ACE_MIND_ID')
batches_last=merge(last,batches_anno,by='ACE_MIND_ID')
cleaned_E=cleaned_E[,batches_last$TechRep]
features=colnames(batches_last)[5:37]
pval_thresh=0.001
lfc_thresh=2

for(feature in features)
{
  batch_subset=batches_last[,c('TechRep','Condition',feature)]
  colnames(batch_subset)=c("TechRep","Condition","feature")
  mod=model.matrix(~feature,data=batch_subset)
  fit=lmFit(2^cleaned_E,mod)
  e=eBayes(fit)
  de=topTable(e,coef=2,sort.by="p",number=nrow(e))
  de$abst=abs(de$t)
  de$feature=row.names(de)
  
  de_sig=de[(abs(de$logFC)>lfc_thresh) & (de$adj.P.Val < pval_thresh),]
  de_sig=de_sig[order(de_sig$abst,decreasing = T),]
  labels=row.names(de_sig)[1:20]
  num_up=sum(as.integer(de_sig$logFC>0))
  num_down=sum(as.integer(de_sig$logFC<0))
  title=paste(feature, "\nup:", num_up,"\ndown:",num_down)
  #make a volcano plot 
  de$padj=de$adj.P.Val
  png(paste(feature,'.png',sep=''),width=7,height=7,units='in',res=300)
  print(volcano_wrapper(de,genes_to_highlight=labels,title=title,pval_thresh=pval_thresh,lfc_thresh=lfc_thresh))
  dev.off()
  #run gProfiler 
  #save sig DE gene list 
  de$feature=feature
  if(dim(de_sig)[1]>0)
  {
    de_sig$feature=feature 
    write.table(de_sig,file=paste(feature,'.de.sig.tsv',sep=''),sep='\t',col.names=T,row.names=T)
    write.table(de,file=paste(feature,'.de.tsv',sep=''),sep='\t',col.names=T,row.names=T)
    foreground=rownames(de_sig)
    background=rownames(de)
  
    out=gProfileR::gprofiler(query = foreground,
                           custom_bg = background,
                           exclude_iea = T,
                           src_filter = c("GO:BP"),
                           organism = "hsapiens",
                           hier_filtering = "moderate")
    out=out[out$term.size<500,]
    
  if(dim(out)[1] >0)
  {  out$feature=feature

    out=out[order(out$p.value,decreasing = T),]
  out$query.number=NULL
  out$significant=NULL
  out$subgraph.number=NULL
  write.table(out,file=paste(feature,'.gProfiler.tsv',sep=''),sep='\t',col.names=T,row.names=F)
  }
  }
}
