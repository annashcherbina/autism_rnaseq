#load the significant genes & clinical feature values 
rm(list=ls())
library(ggplot2)
data=read.table("../../post_processing/NPC_only.corrected_tpm.txt",header=TRUE,check.names = F)

batches=read.table("../../post_processing/merged_rsem/batches.txt",header=TRUE,check.names=F)
batches=batches[batches$Cell=="NPC",]
batches$Sample=factor(batches$Sample)

first=as.data.frame(scale(read.table('../first_vals.txt',header=T,sep='\t',row.names=1,check.names = F)))
first$ACE_MIND_ID=row.names(first)
last=as.data.frame(scale(read.table('../last_vals.txt',header=T,sep='\t',row.names=1,check.names = F)))
last$ACE_MIND_ID=row.names(last)
mapfile=read.table('../APP iPSC pilot study_ID numbers.tsv',header=T,sep='\t')
batches_anno=merge(batches,mapfile,by.y ='iPSC_ID',by.x='Sample')
batches_first=merge(first,batches_anno,by='ACE_MIND_ID')
batches_last=merge(last,batches_anno,by='ACE_MIND_ID')
features=colnames(batches_last)[2:37]

#get the significantly differential genes 
corvals_last=list()
for(comparison in c('ConditionASDDM_vs_ConditionASDN','ConditionASDDM_vs_ConditionTDN','ConditionASDN_vs_ConditionTDN'))
{
print(comparison)
corvals_last[[comparison]]=list()
de=rownames(read.table(paste0('../../post_processing/diff_genes/',comparison,'.de.sig.varfilter.tsv'),header=T,sep='\t',row.names=1))
data_de=data[rownames(data) %in% de,]
gene_names_de=rownames(data_de)
data_de=data_de[,batches_last$TechRep]
for(feature in features)
{
  print(feature)
  corvals_last[[comparison]][[feature]]=list()
  for(gene in gene_names_de)
  {
    corval=cor.test(as.numeric(data_de[gene,]),batches_last[[feature]], method='spearman')
    if(corval$p.value < 0.05)
    {
      corvals_last[[comparison]][[feature]][[gene]]=as.numeric(corval$estimate)
    }
    else
    {
      corvals_last[[comparison]][[feature]][[gene]]=0
      
    }
  }
}
}

#get corvals for first measurement
corvals_first=list()
for(comparison in c('ConditionASDDM_vs_ConditionASDN','ConditionASDDM_vs_ConditionTDN','ConditionASDN_vs_ConditionTDN'))
{
  print(comparison)
  corvals_first[[comparison]]=list()
  de=rownames(read.table(paste0('../../post_processing/diff_genes/',comparison,'.de.sig.varfilter.tsv'),row.names=1,header=T,sep='\t'))
  data_de=data[rownames(data) %in% de,]
  gene_names_de=rownames(data_de)
  data_de=data_de[,batches_first$TechRep]
  for(feature in features)
  {
    print(feature)
    corvals_first[[comparison]][[feature]]=list()
    for(gene in gene_names_de)
    {
      corval=cor.test(as.numeric(data_de[gene,]),batches_first[[feature]],method='spearman')
      if(corval$p.value < 0.05)
      {
        corvals_first[[comparison]][[feature]][[gene]]=as.numeric(corval$estimate)
      }
      else
      {
        corvals_first[[comparison]][[feature]][[gene]]=0
      }
    }
  }
}

#get data frames of feature x gene 
library(ComplexHeatmap)
library(ggplot2)
all_corvals=list() 
all_corvals[['first']]=corvals_first
all_corvals[['last']]=corvals_last
library(pheatmap)
library(matrixStats)
for(comparison in c('ConditionASDDM_vs_ConditionASDN','ConditionASDDM_vs_ConditionTDN','ConditionASDN_vs_ConditionTDN'))
{
  #for(measurement in c('first','last'))
  for(measurement in c('last'))
  {
    cur_corvals=as.data.frame(t(do.call(rbind.data.frame, all_corvals[[measurement]][[comparison]])))
    cur_corvals=cur_corvals[,colMaxs(abs(as.matrix(cur_corvals)))>0.7]
    cur_corvals=cur_corvals[rowMaxs(abs(as.matrix(cur_corvals)))>0.7,]
    obs_string=paste(measurement,'observation')
    print(dim(cur_corvals))
    pheatmap(as.matrix(cur_corvals),main=paste('Spearman rho', comparison, obs_string,sep='\n'), fontsize_row=2)
    #png(paste0('spearman.',comparison,'.',measurement,'.png'),width=30,height=6,res=160,units='in')
    #print(Heatmap(as.matrix(cur_corvals),name=paste('Spearman cor', comparison, obs_string,sep='\n'),
    #        column_title='Gene',
    #        row_title = 'Clinical Feature',
    #        column_title_side = 'bottom',
    #        row_title_side = 'right'))
    #dev.off()
    
    write.table(cur_corvals,file=paste0('spearman.',comparison,'.',measurement,'.tsv'),
                row.names=T,
                col.names=T,
                sep='\t')
  }
}

