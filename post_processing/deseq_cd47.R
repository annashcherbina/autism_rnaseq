rm(list=ls())
library(limma)
library(sva)
library(ggplot2)
library(DESeq2)
library("BiocParallel")
parallelFlag=TRUE

data=read.table("merged_rsem/rna.counts.txt.expected_count",header=TRUE)
batches=read.table("merged_rsem/batches.txt",header=TRUE)

# columns to paste together
cols <- c( 'GeneID' , 'GeneName')

# create a new column `x` with the three columns collapsed together
rownames(data) <- apply( data[ , cols ] , 1 , paste , collapse = "-" )

# remove the unnecessary columns
data <- data[ , !( names( data ) %in% cols ) ]

#remove genes with 0 counts 
data=data[rowSums(data)>0,]
#DESIGN MATRIX
mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~Condition+Cell,data=batches)

#RUN SVA 
sva.obj=svaseq(as.matrix(data),mod1,mod0)
sur_var=data.frame(sva.obj$sv)
batches=cbind(batches,sur_var)

#CREATE A NEW GROUP THAT ACCOUNTS FOR DEL/DUP/WT + MACROCEPHALY 
batches$Group=paste(batches$Cell,batches$Condition,sep='_')

#Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(data),
                              colData = batches,
                              design = ~Group+X1+X2+X3+X4+X5+X6+X7+X8+X9)

#Run the differential analysis
dds <- DESeq(dds,parallel = TRUE)
res=results(dds)
summary(res)

res=results(dds,independentFiltering=FALSE)
summary(res)

res=results(dds,filterFun = ihw)
summary(res)

group1=c("IPSC_ASDDM",
         "IPSC_ASDDM",
         "IPSC_ASDN",
         "NPC_ASDDM",
         "NPC_ASDDM",
         "NPC_ASDN")
group2=c("IPSC_ASDN",
         "IPSC_TDN",
         "IPSC_TDN",
         "NPC_ASDN",
         "NPC_TDN",
         "NPC_TDN")

comparisons=c("IPSC_ASDDM_ASDN",
              "IPSC_ASDDM_TDN",
              "IPSC_ASDN_TDN",
              "NPC_ASDDM_ASDN",
              "NPC_ASDDM_TDN",
              "NPC_ASDN_TDN")

##get the results for the various contrasts 
numcomparisons=length(comparisons)
for(i in seq(1,numcomparisons))
{
  res=results(dds, contrast=c("Group", group1[i],group2[i]),parallel=TRUE)
  res$logPadj=-1*log10(res$padj)
  res=as.data.frame(res)
  res=na.omit(res)
  print(comparisons[i])
  print(res['ENSG00000196776.16-CD47',])
}