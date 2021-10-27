rm(list=ls())
library(limma)
library(sva)
library(ggplot2)

data=read.table("merged_rsem/rna.counts.txt.tpm",header=TRUE)
batches=read.table("merged_rsem/batches.txt",header=TRUE)
batches$Sample=factor(batches$Sample)

# columns to paste together
cols <- c( 'GeneID' , 'GeneName')

# create a new column `x` with the three columns collapsed together
rownames(data) <- apply( data[ , cols ] , 1 , paste , collapse = "||" )

# remove the unnecessary columns
data <- data[ , !( names( data ) %in% cols ) ]

#WE RESTRICT OUT ANALLYSIS TO NPC, REMOVE ANY IPSC DATA COLUMNS 
batches=batches[batches$Cell=='NPC',]
data=subset(data,select=c(batches$TechRep))

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
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: PC1 vs PC2")+
  theme_bw()

ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,shape=pca_df$Condition,color=pca_df$Sample))+
  geom_point(size=3)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: PC1 vs PC2")+
  theme_bw()



ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,color=pca_df$Batch))+
  geom_point(size=3)+
  xlab(paste("PC2:",var_explained[2]))+
  ylab(paste("PC3:",var_explained[3]))+
  ggtitle("PCA: PC2 vs PC3")+
  theme_bw()

#DESIGN MATRIX
mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~Condition,data=batches)

#RUN SVA 
sva.obj=sva(E,mod1,mod0)
sur_var=data.frame(sva.obj$sv)
cleaned_E=removeBatchEffect(E,covariates=sur_var, batch2=batches$Batch,design=mod1)
data.pca=prcomp(t(cleaned_E))
var_explained=as.character(round(100*data.pca$sdev^2/sum(data.pca$sdev^2),2))
barplot(100*data.pca$sdev^2/sum(data.pca$sdev^2),las=2,ylab="% Variance Explained",xlab="Principal Component")

pca_df=data.frame(data.pca$x)
pca_df=cbind(pca_df,batches)

ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,color=pca_df$Batch,label=pca_df$TechRep))+
  geom_point(size=3)+
  geom_text(hjust=0, vjust=0, size=3)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: SVA correction")+
  theme_bw()

ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,color=pca_df$Sample))+
  geom_point(size=3)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: SVA correction")+
  theme_bw()

ggplot(data=pca_df,aes(x=pca_df$PC1,y=pca_df$PC2,color=pca_df$Condition))+
  geom_point(size=3)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: SVA correction")+
  theme_bw()

#LIMMA REGRESSION ANALYSIS
mod2=model.matrix(~0+Condition,data=batches)
fit=lmFit(2^cleaned_E,mod2)
cont.matrix=makeContrasts(ASDDM_ASDN="ConditionASDDM-ConditionASDN",
                          ASDDM_TDN="ConditionASDDM-ConditionTDN",
                          ASDN_TDN="ConditionASDN-ConditionTDN",levels=mod2)

fit2=contrasts.fit(fit,cont.matrix)
e=eBayes(fit2)
comparisons=colnames(cont.matrix)
for(i in seq(1,length(comparisons)))
{
  print(comparisons[i])
  write.table(topTable(e,number=nrow(e),coef=i),paste(file=comparisons[i],'diff','tsv',sep='.'),row.names=TRUE,col.names=TRUE,sep='\t')
  volcanoplot(e,coef=i,style='p-value',highlight=10)
}

tpm_cd99=tpm['ENSG00000223773.7||CD99P1',]
batches$cd99=tpm_cd99
ggplot(batches,
       aes(x=batches$Condition,
           y=batches$cd99))+
  geom_boxplot()+
  xlab("Condition")+
  ylab("Corrected TPM cd99")
write.table(batches,file="tpm_cd47.txt",row.names = TRUE,col.names = TRUE,sep='\t')

tpm=2^cleaned_E
write.table(tpm,file="NPC_only.corrected_tpm.txt",row.names=TRUE,col.names=TRUE,sep='\t')
tpm_cd47=tpm['ENSG00000196776.16||CD47',]
batches$cd47=tpm_cd47
ggplot(batches,
       aes(x=batches$Condition,
           y=batches$cd47))+
  geom_boxplot()+
  xlab("Condition")+
  ylab("Corrected TPM CD47")
write.table(batches,file="tpm_cd47.txt",row.names = TRUE,col.names = TRUE,sep='\t')

tpm_pax6=tpm['ENSG00000007372.23||PAX6',]
batches$pax6=tpm_pax6
ggplot(batches,
       aes(x=batches$Condition,
           y=batches$pax6))+
  geom_boxplot()+
  xlab("Condition")+
  ylab("Corrected TPM pax6")

write.table(batches,file="tpm_pax6.txt",row.names = TRUE,col.names = TRUE,sep='\t')
