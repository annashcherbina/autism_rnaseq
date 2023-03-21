rm(list=ls())
library(ggplot2)

genes_of_interest=c("CD99",
                    "IER3")
data=read.table("NPC_only.corrected_tpm.txt",header=TRUE,sep='\t',row.names = 1,check.names=FALSE)

#subset data to genes of interest
data=data[genes_of_interest,]
data=as.data.frame(t(data))
batches=read.table("merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)
data=merge(data,batches,by=0)
data$Replicate=data$Row.names
data$Row.names=NULL
data=data[order(data$CD99),]
data$Replicate=factor(data$Replicate,levels=data$Replicate)
p1=ggplot(data=data,
       aes(x=data$Replicate,
           y=data$CD99))+geom_bar(stat='identity')+xlab("Replicate")+ylab("TPM CD99")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2=ggplot(data=data,
          aes(x=data$Replicate,
              y=data$IER3))+geom_bar(stat='identity')+xlab("Replicate")+ylab("TPM IER3")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

source("~/helpers.R")
multiplot(p1,p2)

