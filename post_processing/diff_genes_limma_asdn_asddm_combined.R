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

#DESIGN MATRIX
mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~Cell+Condition,data=batches)

#RUN SVA 
sva.obj=sva(E,mod1,mod0)
sur_var=data.frame(sva.obj$sv)
cleaned_E=removeBatchEffect(E,covariates=sur_var, batch2=batches$Pool,design=mod1)

#CREATE A NEW GROUP THAT ACCOUNTS FOR DEL/DUP/WT + MACROCEPHALY 
batches$ASD="ASD"
batches$ASD[batches$Condition=="TDN"]="TDN"
batches$Group=paste(batches$Cell,batches$ASD,sep='_')
#LIMMA REGRESSION ANALYSIS
mod2=model.matrix(~0+Group,data=batches)
fit=lmFit(2^cleaned_E,mod2)

cont.matrix=makeContrasts(IPSC="GroupIPSC_ASD-GroupIPSC_TDN",
                          NPC="GroupNPC_ASD-GroupNPC_TDN",levels=mod2)

fit2=contrasts.fit(fit,cont.matrix)
e=eBayes(fit2)
comparisons=colnames(cont.matrix)
for(i in seq(1,length(comparisons)))
{
  print(comparisons[i])
  write.table(topTable(e,number=nrow(e),coef=i),paste(file=comparisons[i],'diff','tsv',sep='.'),row.names=TRUE,col.names=TRUE,sep='\t')
  volcanoplot(e,coef=i,style='p-value',highlight=10)
}
