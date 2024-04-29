rm(list=ls())
library(ggplot2)

genes_of_interest=c("HSPB1")
data=read.table("merged_rsem/tpm.txt",header=TRUE,sep='\t',row.names = 1,check.names=FALSE)

#subset data to genes of interest
data=data[data$Symbol %in% c(genes_of_interest),]
data$gene_id=NULL
data$gene_symbol=NULL
data$Symbol=NULL
data=colSums(data)
batches=read.table("merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)
data=merge(data,batches,by=0)
data$Replicate=data$Row.names
data$Row.names=NULL
data$HSP27=data$x
data=data[order(data$HSP27),]
data$Replicate=factor(data$Replicate,levels=data$Replicate)
rownames(data)=data$Replicate
p1=ggplot(data=data,
       aes(x=data$Replicate,
           y=data$HSP27))+geom_bar(stat='identity')+xlab("Replicate")+ylab("TPM HSP27")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

source("~/helpers.R")
#multiplot(p1,p2)

