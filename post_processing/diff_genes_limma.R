rm(list=ls())
library(limma)
library(sva)
library(ggplot2)

data=read.table("merged_rsem/rna.counts.txt.tpm",header=TRUE)
batches=read.table("merged_rsem/batches.txt",header=TRUE)

# columns to paste together
cols <- c( 'GeneID' , 'GeneName')

# create a new column `x` with the three columns collapsed together
rownames(data) <- apply( data[ , cols ] , 1 , paste , collapse = "-" )

# remove the unnecessary columns
data <- data[ , !( names( data ) %in% cols ) ]

#remove genes with 0 counts 
data=data[rowSums(data)>0,]
#get TPM, quantile normalize
tpm=voom(data,normalize.method = 'quantile')
E=tpm$E
E=round(E,2)

data.pca=prcomp(t(E))
var_explained=as.character(round(100*data.pca$sdev^2/sum(data.pca$sdev^2),2))
barplot(100*data.pca$sdev^2/sum(data.pca$sdev^2),las=2,ylab="% Variance Explained",xlab="Principal Component")

pca_df=data.frame(data.pca$x)
pca_df=cbind(pca_df,batches)
pca_df$Sample=factor(pca_df$Sample)
ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,shape=pca_df$Cell,color=pca_df$Batch))+
  geom_point(size=3)+
  xlab("PC1: 57.91%")+
  ylab("PC2: 4.92%")+
  ggtitle("PCA: PC1 vs PC2")+
  theme_bw()

ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,shape=pca_df$Cell,color=pca_df$Sample))+
  geom_point(size=3)+
  xlab("PC1: 57.91%")+
  ylab("PC2: 4.92%")+
  ggtitle("PCA: PC1 vs PC2")+
  theme_bw()


#DESIGN MATRIX
mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~Cell+Condition,data=batches)

#RUN SVA 
sva.obj=sva(E,mod1,mod0)
sur_var=data.frame(sva.obj$sv)
cleaned_E=removeBatchEffect(E,covariates=sur_var, batch2=batches$Pool,design=mod1)
data.pca=prcomp(t(cleaned_E))
var_explained=as.character(round(100*data.pca$sdev^2/sum(data.pca$sdev^2),2))
barplot(100*data.pca$sdev^2/sum(data.pca$sdev^2),las=2,ylab="% Variance Explained",xlab="Principal Component")

pca_df=data.frame(data.pca$x)
pca_df=cbind(pca_df,batches)


ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,shape=pca_df$Cell,color=pca_df$Batch,label=pca_df$TechRep))+
  geom_point(size=3)+
  geom_text(hjust=0, vjust=0, size=3)+
  xlab("PC1: 70.03%")+
  ylab("PC2: 1.8%")+
  ggtitle("PCA: SVA correction")+
  theme_bw()

ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,shape=pca_df$Cell,color=pca_df$Condition,label=pca_df$TechRep))+
  geom_point(size=3)+
  geom_text(hjust=0, vjust=0, size=3)+
  xlab("PC1: 70.03%")+
  ylab("PC2: 1.8%")+
  ggtitle("PCA: SVA correction")+
  theme_bw()

ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC3,shape=pca_df$Cell,color=pca_df$Condition,label=pca_df$TechRep))+
  geom_point(size=3)+
  geom_text(hjust=0, vjust=0, size=3)+
  xlab("PC1: 70.03%")+
  ylab("PC3: 1.52%")+
  ggtitle("PCA: SVA correction")+
  theme_bw()

ggplot(data=pca_df,aes(x=pca_df$PC2,y=pca_df$PC3,shape=pca_df$Cell,color=pca_df$Condition,label=pca_df$TechRep))+
  geom_point(size=3)+
  geom_text(hjust=0, vjust=0, size=3)+
  xlab("PC2: 1.8%")+
  ylab("PC3: 1.52%")+
  ggtitle("PCA: SVA correction")+
  theme_bw()


#CREATE A NEW GROUP THAT ACCOUNTS FOR DEL/DUP/WT + MACROCEPHALY 
batches$Group=paste(batches$Cell,batches$Condition,sep='_')
#LIMMA REGRESSION ANALYSIS
mod2=model.matrix(~0+Group,data=batches)
fit=lmFit(2^cleaned_E,mod2)
cont.matrix=makeContrasts(IPSC_ASDDM_ASDN="GroupIPSC_ASDDM-GroupIPSC_ASDN",
                          IPSC_ASDDM_TDN="GroupIPSC_ASDDM-GroupIPSC_TDN",
                          IPSC_ASDN_TDN="GroupIPSC_ASDN-GroupIPSC_TDN",
                          NPC_ASDDM_ASDN="GroupNPC_ASDDM-GroupNPC_ASDN",
                          NPC_ASDDM_TDN="GroupNPC_ASDDM-GroupNPC_TDN",
                          NPC_ASDN_TDN="GroupNPC_ASDN-GroupNPC_TDN",levels=mod2)

fit2=contrasts.fit(fit,cont.matrix)
e=eBayes(fit2)
comparisons=colnames(cont.matrix)
for(i in seq(1,length(comparisons)))
{
  print(comparisons[i])
  write.table(topTable(e,number=nrow(e),coef=i),paste(file=comparisons[i],'diff','tsv',sep='.'),row.names=TRUE,col.names=TRUE,sep='\t')
  volcanoplot(e,coef=i,style='p-value',highlight=10)
}

tpm=2^cleaned_E
write.table(tpm,file="corrected_tpm.txt",row.names=TRUE,col.names=TRUE,sep='\t')
tpm_cd47=tpm['ENSG00000196776.16-CD47',]
batches$cd47=tpm_cd47
ggplot(batches,
       aes(x=batches$Condition,
           y=batches$cd47,
           fill=batches$Cell))+
  geom_boxplot()
write.table(batches,file="tpm_cd47.txt",row.names = TRUE,col.names = TRUE,sep='\t')

tpm_pax6=tpm['ENSG00000007372.23-PAX6',]
batches$pax6=tpm_pax6
write.table(batches,file="tpm_pax6.txt",row.names = TRUE,col.names = TRUE,sep='\t')
