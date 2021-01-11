library(ggplot2)
source("~/helpers.R")
data=t(read.table("genes_of_interest.tpm.txt",header=TRUE,sep='\t',row.names = 1,check.names=FALSE))
batches=read.table("merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)
data=merge(data,batches,by=0)
data$Condition=factor(data$Condition,levels=c("TDN","ASDN","ASDDM"))
data1=data[data$Cell=='IPSC',]
data2=data[data$Cell=="NPC",]
p1=ggplot(data=data1,
       aes(x=Condition,y=`ENSG00000007372.23-PAX6`))+
    geom_boxplot(color='#FF0000',outlier.shape = NA)+
    geom_jitter()+
    theme_bw(15)+
  ylab("TPM")+
  ggtitle("PAX6 TPM, IPSC")
p2=ggplot(data=data2,
          aes(x=Condition,y=`ENSG00000007372.23-PAX6`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("PAX6 TPM, NPC")

p3=ggplot(data=data1,
          aes(x=Condition,y=`ENSG00000106689.11-LHX2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("LHX2 TPM, IPSC")
p4=ggplot(data=data2,
          aes(x=Condition,y=`ENSG00000106689.11-LHX2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("LHX2 TPM, NPC")

p5=ggplot(data=data1,
          aes(x=Condition,y=`ENSG00000138083.5-SIX3`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("SIX3 TPM, IPSC")
p6=ggplot(data=data2,
          aes(x=Condition,y=`ENSG00000138083.5-SIX3`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("SIX3 TPM, NPC")

p7=ggplot(data=data1,
          aes(x=Condition,y=`ENSG00000170370.12-EMX2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("EMX2 TPM, IPSC")
p8=ggplot(data=data2,
          aes(x=Condition,y=`ENSG00000170370.12-EMX2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("EMX2 TPM, NPC")
p9=ggplot(data=data1,
          aes(x=Condition,y=`ENSG00000106689.11-LHX2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("LHX2 TPM, IPSC")
p10=ggplot(data=data2,
          aes(x=Condition,y=`ENSG00000181449.4-SOX2`))+
  geom_boxplot(color='#FF0000',outlier.shape = NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("SOX2 TPM, NPC")

multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,cols=5)
