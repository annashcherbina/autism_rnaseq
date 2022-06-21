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

#DESIGN MATRIX
mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~Condition,data=batches)

#RUN SVA 
sva.obj=sva(E,mod1,mod0)
sur_var=data.frame(sva.obj$sv)
cleaned_E=removeBatchEffect(E,covariates=sur_var, batch2=batches$Batch,design=mod1)

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
  batch_subset$feature=sample(batches$feature,length(batches$feature),replace=FALSE)
  colnames(batch_subset)=c("TechRep","Condition","feature")
  mod=model.matrix(~feature,data=batch_subset)
  fit=lmFit(2^cleaned_E,mod)
  e=eBayes(fit)
  de=topTable(e,coef=2,sort.by="p",number=nrow(e))
  de$feature=row.names(de)
  de$padj=de$adj.P.Val
  de$feature=feature
  write.table(de,file=paste(feature,'.de.permuted.tsv',sep=''),sep='\t',col.names=T,row.names=T)
}
