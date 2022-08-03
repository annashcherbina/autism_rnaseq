library(ggplot2)
source("~/helpers.R")
data=t(read.table("genes_of_interest_simon.tpm.txt",header=TRUE,sep='\t',row.names = 1,check.names=FALSE))
genes=colnames(data)
batches=read.table("../merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)
batches=batches[batches$Cell=="NPC",]

data=merge(data,batches,by=0)
data$Condition=factor(data$Condition,levels=c("TDN","ASDN","ASDDM"))
plots=list()
for(gene in genes)
{
  plots[[gene]]=ggplot(data=data,
                       aes(x=Condition,y=data[[gene]]))+
    geom_boxplot(color='#FF0000',outlier.shape = NA)+
    geom_jitter()+
    theme_bw(8)+
    ylab("TPM")+
    ggtitle(paste(gene,"TPM, NPC"))
  
}
png('genes_of_interest.png',height=10,width=8,units='in',res=160)
multiplot(plotlist = plots, cols = 5)
dev.off()

