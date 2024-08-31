rm(list=ls())
library(ggplot2)

#genes_of_interest=read.table('cell_stress.txt')$V1
#genes_of_interest=read.table('heat_shock.txt')$V1
#genes_of_interest=read.table('iron_proteins.txt')$V1
#genes_of_interest=read.table('selenocysteine.txt')$V1
#genes_of_interest=read.table('iron_proteins.txt')$V1
#genes_of_interest=read.table('cellular_lipid_metabolic_process.txt')$V1
#genes_of_interest=read.table('oxidation.txt')$V1
genes_of_interest=c("PKM")
data=read.table("../NPC_only.corrected_tpm.txt",header=TRUE,sep='\t',row.names = 1,check.names=FALSE)
genes_of_interest=intersect(rownames(data),genes_of_interest)
#subset data to genes of interest
data=data[genes_of_interest,]
data$Symbol=NULL
data$gene_id=NULL
data=as.data.frame(t(data))
batches=read.table("../merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)

data=merge(data,batches,by=0)

conditions = c("TDN","ASDN","ASDDM")
data=data[data$Cell=='NPC',]
data$Condition=factor(data$Condition, levels=conditions)
data$Condition=as.character(data$Condition)
data$Condition[data$Condition=='TDN']='TD-N'
data$Condition[data$Condition=='ASDDM']='ASD-DM'
data$Condition[data$Condition=='ASDN']='ASD-N'

plot.list = lapply(genes_of_interest, function(gene){
    p <- ggplot(data=data,
                aes(x=data$Condition,
                    y=data[[gene]]))+
      geom_boxplot(color='#FF0000',outlier.shape=NA)+
      geom_jitter()+
      ylab("TPM")+
      xlab("Condition")+
      ggtitle(paste(gene,'\nTPM,',"NPC"))+
      theme_bw(10)
    return(p)
  })
library(gridExtra)
margin = theme(plot.margin = unit(c(0.25,0.25,1,1), "cm"))
#grid.arrange(grobs = lapply(plot.list, ncol=3,nrow=3)
grid.arrange(grobs = lapply(plot.list, "+", margin),ncol=1,nrow=1)
