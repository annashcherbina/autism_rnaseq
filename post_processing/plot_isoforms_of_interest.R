rm(list=ls())
library(ggplot2)
library(reshape2)
source("~/helpers.R")
genes_of_interest=c("ISG20",
                    "APOBEC3H",
                    "IFIH1",
                    "IFIT3",
                    "CLDN6",
                    "IFIT1",
                    "KLF11",
                    "MX1",
                    "FOXJ1",
                    "HERC5",
                    "PARP14",
                    "SPINK5",
                    "LUM",
                    "UCP3",
                    "RHOU",
                    "WNT7B")
genes_of_interest=c("IFI6",
                    "HTRA1",
                    "IER3",
                    "TFPI2", 
                    "GDF15", 
                    "GAPDH",
                    "HSP90AB1",
                    "RPL30", 
                    "RPS17", 
                    "PPIA", 
                    "POP4",
                    "KCTD13", 
                    "ALDOA",
                    "MAPK3",
                    "TAOK2")
genes_of_interest=c("ANGPT2","DOC2A","ISM2","TBX6")
genes_of_interest=c("IER3","PKIB","PKM","CLU")
genes_of_interest=c("GRN")
data=read.table("merged_rsem/tpm.isoform.txt",header=TRUE,sep='\t',row.names = 1,check.names=F)
pct=read.table("merged_rsem/iso_pct.isoform.txt",header=T,sep='\t',row.names=1,check.names=F)
batches=read.table("merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)
genes=data$gene_id
symbols=data$Symbol
batches=batches[batches$Cell=="NPC",]
data=data[,rownames(batches)]
data$Symbol=symbols 
pct=pct[,rownames(batches)]
pct$Symbol=symbols 

#subset data to genes of interest
data=data[data$Symbol %in% genes_of_interest,]
pct=pct[pct$Symbol %in% genes_of_interest,]

tpm_plots=list()
iso_pct_plots=list()
for(gene in genes_of_interest)
{
  data_sub=data[data$Symbol==gene,]
  pct_sub=pct[pct$Symbol==gene,]
  data_sub$Symbol=NULL
  pct_sub$Symbol=NULL
  data_sub$isoform=rownames(data_sub)
  pct_sub$isoform=rownames(pct_sub)
  data_sub=melt(data_sub)
  pct_sub=melt(pct_sub)
  
  data_sub$Condition=batches[data_sub$variable,'Condition']
  pct_sub$Condition=batches[pct_sub$variable,'Condition']
  
  p1=ggplot(data_sub,
         aes(x=Condition,
             y=value,
             fill=isoform))+geom_boxplot()+ylab("TPM")+ggtitle(gene)
  p2=ggplot(pct_sub,
            aes(x=Condition,
                y=value,
                fill=isoform))+geom_boxplot()+ylab("% of gene reads")+ggtitle(gene)
  png(paste0(gene,".png"))
  multiplot(p1,p2,cols=1)
  dev.off()
    
}