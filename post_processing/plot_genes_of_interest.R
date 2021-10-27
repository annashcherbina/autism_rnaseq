library(ggplot2)
source("~/helpers.R")
data=t(read.table("genes_of_interest.tpm.txt",header=TRUE,sep='\t',row.names = 1,check.names=FALSE))
batches=read.table("merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)
batches=batches[batches$Cell=="NPC",]

data=merge(data,batches,by=0)
data$Condition=factor(data$Condition,levels=c("TDN","ASDN","ASDDM"))
p1=ggplot(data=data,
       aes(x=Condition,y=`ENSG00000007372.23-PAX6`))+
    geom_boxplot(color='#FF0000',outlier.shape = NA)+
    geom_jitter()+
    theme_bw(15)+
  ylab("TPM")+
  ggtitle("PAX6 TPM, NPC")

p2=ggplot(data=data,
          aes(x=Condition,y=`ENSG00000106689.11-LHX2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("LHX2 TPM, NPC")

p3=ggplot(data=data,
          aes(x=Condition,y=`ENSG00000138083.5-SIX3`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("SIX3 TPM, NPC")

p4=ggplot(data=data,
          aes(x=Condition,y=`ENSG00000170370.12-EMX2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("EMX2 TPM, NPC")
p9=ggplot(data=data,
          aes(x=Condition,y=`ENSG00000106689.11-LHX2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("LHX2 TPM, NPC")

p5=ggplot(data=data,
          aes(x=Condition,y=`ENSG00000196776.16-CD47`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("CD47 TPM, NPC")


p6=ggplot(data=data,
          aes(x=Condition,y=`ENSG00000002586.20-CD99`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("CD99 TPM, NPC")

p7=ggplot(data=data,
          aes(x=Condition,y=`ENSG00000102181.21-CD99L2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("CD99L2 TPM, NPC")


p8=ggplot(data=data,
          aes(x=Condition,y=`ENSG00000223773.7-CD99P1`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("CD99P1 TPM, NPC")


multiplot(p1,p2,p3,p4,p5,p6,p7,p8,cols=4)
