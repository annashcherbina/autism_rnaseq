#load the significant genes & clinical feature values 
rm(list=ls())
library(ggplot2)
data=read.table("../post_processing/corrected_tpm.txt",header=TRUE,check.names = F)
gene_names=as.data.frame(str_split_fixed(rownames(data),'-',2))
data$Gene=gene_names$V2

batches=read.table("../post_processing/merged_rsem/batches.txt",header=TRUE,check.names=F)
batches$Sample=factor(batches$Sample)

first=as.data.frame(scale(read.table('first_vals.txt',header=T,sep='\t',row.names=1,check.names = F)))
first$ACE_MIND_ID=row.names(first)
last=as.data.frame(scale(read.table('last_vals.txt',header=T,sep='\t',row.names=1,check.names = F)))
last$ACE_MIND_ID=row.names(last)
mapfile=read.table('APP iPSC pilot study_ID numbers.tsv',header=T,sep='\t')
batches_anno=merge(batches,mapfile,by.y ='iPSC_ID',by.x='Sample')
batches_first=merge(first,batches_anno,by='ACE_MIND_ID')
batches_last=merge(last,batches_anno,by='ACE_MIND_ID')
features=colnames(batches_last)[2:37]

#get the significantly differential genes 
corvals_last=list()
for(comparison in c('ASDDM_ASDN','ASDDM_TDN','ASDN_TDN'))
{
print(comparison)
corvals_last[[comparison]]=list()
de=read.table(paste0('all.sig.',comparison,'.tsv'),header=F,sep='\t')$V1
data_de=data[data$Gene %in% de,]
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
for(comparison in c('ASDDM_ASDN','ASDDM_TDN','ASDN_TDN'))
{
  print(comparison)
  corvals_first[[comparison]]=list()
  de=read.table(paste0('all.sig.',comparison,'.tsv'),header=F,sep='\t')$V1
  data_de=data[data$Gene %in% de,]
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
for(comparison in c('ASDDM_ASDN','ASDDM_TDN','ASDN_TDN'))
{
  for(measurement in c('first','last'))
  {
    cur_corvals=as.data.frame(t(do.call(rbind.data.frame, all_corvals[[measurement]][[comparison]])))
    cur_corvals=cur_corvals[,colSums(cur_corvals)>0.05]
    cur_corvals=cur_corvals[rowSums(cur_corvals)>0.05,]
    gene_names=as.data.frame(str_split_fixed(rownames(cur_corvals),'\\.',3))
    rownames(cur_corvals)=gene_names$V3
    obs_string=paste(measurement,'observation')
    print(dim(cur_corvals))
    png(paste0('spearman.',comparison,'.',measurement,'.png'),width=20,height=6,res=160,units='in')
    print(Heatmap(t(as.matrix(cur_corvals)),name=paste('Spearman cor', comparison, obs_string,sep='\n'),
            column_title='Gene',
            row_title = 'Clinical Feature',
            column_title_side = 'bottom',
            row_title_side = 'right'))
    dev.off()
    
    write.table(cur_corvals,file=paste0('spearman.',comparison,'.',measurement,'.tsv'),
                row.names=T,
                col.names=T,
                sep='\t')
  }
}

