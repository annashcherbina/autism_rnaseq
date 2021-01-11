library(ggplot2)
data=read.table("tpm_pax6.txt.txt",header=TRUE,sep='\t')
data1=data[data$Cell=='IPSC',]
data2=data[data$Cell=="NPC",]
p1=ggplot(data=data1,
       aes(x=Condition,y=cd47))+
    geom_boxplot()+
    geom_jitter()+
    theme_bw(15)+
  ylab("TPM")+
  ggtitle("CD47 TPM, IPSC")
p2=ggplot(data=data2,
          aes(x=Condition,y=cd47))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw(15)+
  ylab("TPM")+
  ggtitle("CD47 TPM, NPC")
source("~/helpers.R")
multiplot(p1,p2,cols=1)
