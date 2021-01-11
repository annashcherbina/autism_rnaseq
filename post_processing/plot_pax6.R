library(ggplot2)
data=read.table("tpm_pax6.txt",header=TRUE,sep='\t')
data$Condition=factor(data$Condition,levels=c("TDN","ASDN","ASDDM"))
data1=data[data$Cell=='IPSC',]
data2=data[data$Cell=="NPC",]
p1=ggplot(data=data1,
          aes(x=Condition,y=pax6))+
  geom_boxplot(color='#FF0000',outlier.shape=NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("PAX6 TPM, IPSC")
p2=ggplot(data=data2,
          aes(x=Condition,y=pax6))+
  geom_boxplot(color='#FF0000',outlier.shape=NA)+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("PAX6 TPM, NPC")
source("~/helpers.R")
multiplot(p1,p2,cols=1)
