rm(list=ls())
library(limma)
library(sva)
library(ggplot2)
library(gProfileR)
source('../volcano_wrapper.R')
data=read.table("../post_processing/merged_rsem/rna.counts.txt.tpm",header=TRUE,check.names = F)
batches=read.table("../post_processing/merged_rsem/batches.txt",header=TRUE)
batches$Sample=factor(batches$Sample)

#sum tpm across genes w/ same name
data=aggregate(data[,3:ncol(data)],by=list(GeneID=data$GeneName),FUN=sum)
row.names(data)=data$GeneID
data$GeneID=NULL

#WE RESTRICT OUT ANALLYSIS TO NPC, REMOVE ANY IPSC DATA COLUMNS 
batches=batches[batches$Cell=='NPC',]
data=subset(data,select=c(batches$TechRep))

#remove genes with 0 counts 
data=data[rowSums(data)>0,]

#remove lncRNA, snoRNA, and similar
genes_to_keep=read.table("../gencode_to_category_filtered.txt",header=T)$g_name
data=data[row.names(data) %in% genes_to_keep,]

#get TPM, quantile normalize
tpm=voom(data,normalize.method = 'quantile')
E=tpm$E
E=round(E,2)

mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~Condition,data=batches)

#RUN SVA 
sva.obj=sva(E,mod1,mod0)
sur_var=data.frame(sva.obj$sv)
cleaned_E=removeBatchEffect(E,covariates=sur_var, batch2=batches$Batch,design=mod1)

#LIMMA REGRESSION ANALYSIS
first=read.table('first_vals.txt',header=T,sep='\t',row.names=1)
first$ACE_MIND_ID=row.names(first)
last=read.table('last_vals.txt',header=T,sep='\t',row.names=1)
last$ACE_MIND_ID=row.names(last)
mapfile=read.table('APP iPSC pilot study_ID numbers.tsv',header=T,sep='\t')
batches_anno=merge(batches,mapfile,by.y ='iPSC_ID',by.x='Sample')
batches_first=merge(first,batches_anno,by='ACE_MIND_ID')
batches_last=merge(last,batches_anno,by='ACE_MIND_ID')
cleaned_E=cleaned_E[,batches_last$TechRep]
features=c("iq_viq",
           "iq_nviq",
           "iq_fsiq",
           "vine_communication_ss",
           "vine_abc_ss",
           "mori_total_l3_occipital_l")
pval_thresh=0.01
lfc_thresh=1

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
                  
  foreground=de_sig$feature
  background=de$feature
  
  out=gProfileR::gprofiler(query = foreground,
                           custom_bg = background,
                           exclude_iea = T,
                           src_filter = c("GO:BP"),
                           organism = "hsapiens",
                           hier_filtering = "moderate")
  out=out[order(out$p.value,decreasing = T),]
  out=out[out$term.size<500,]
  out$query.number=NULL
  out$significant=NULL
  out$subgraph.number=NULL
  out$feature=feature
  write.table(out,file=paste(feature,'.gProfiler.tsv',sep=''),sep='\t',col.names=T,row.names=F)
  #save sig DE gene list 
  de$feature=feature
  de_sig$feature=feature 
  write.table(de_sig,file=paste(feature,'.de.sig.tsv',sep=''),sep='\t',col.names=T,row.names=T)
  write.table(de,file=paste(feature,'.de.tsv',sep=''),sep='\t',col.names=T,row.names=T)
}
