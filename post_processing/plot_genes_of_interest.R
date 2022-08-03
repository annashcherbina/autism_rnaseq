rm(list=ls())
library(ggplot2)

genes_of_interest=c("HERC5",
                    "IFIT3",
                    "IFITM2",
                    "ISG15",
                    "ISG20",
                    "MX1",
                    "CD99",
                    "KIF22",
                    "MAZ",
                    "PLCG1",
                    "FLNB",
                    "PPP4C",
                    "DOC2A",
                    "PRRT2",
                    "SEZ6L2")
#genes_of_interest=c("GAPDH","TFRC","HSP90AB1","RPL30","RPS17","RPL37A","PPIA","RNA18S1","B2M","RPLPO","HPRT1","POP4","CDKN1B","ELF1")
genes_of_interest=c("CLDN6","IFIT1","KLF11", "MX1","FOXJ1","HERC5","PARP14","ISG20","SPINK5","LUM","UCP3","ETS1","SSTR2","IFIH1","C5AR1")
data=read.table("NPC_only.corrected_tpm.txt",header=TRUE,sep='\t',row.names = 1,check.names=FALSE)
#subset data to genes of interest
data=data[genes_of_interest,]
data=as.data.frame(t(data))
batches=read.table("merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)

data=merge(data,batches,by=0)

conditions = c("TDN","ASDN","ASDDM")

data$Condition=factor(data$Condition, levels=conditions)


plot.list = lapply(genes_of_interest, function(gene){
    p <- ggplot(data=data,
                aes(x=data$Condition,
                    y=data[[gene]]))+
      geom_boxplot(color='#FF0000',outlier.shape=NA)+
      geom_jitter()+
      ylab("TPM")+
      xlab("Condition")+
      ggtitle(paste(gene,'\nTPM,',"NPC"))+
      theme_bw(20)
    return(p)
  })
library(gridExtra)
grid.arrange(grobs = plot.list, ncol=5,nrow=3)
