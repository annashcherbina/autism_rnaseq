library(DESeq2)
data=read.table("merged_rsem/rna.counts.txt.expected_count",header=TRUE)
data$GeneName=NULL

batches=read.table("merged_rsem/batches.txt",header=TRUE)
batches$Sample=factor(batches$Sample)

#WE RESTRICT OUT ANALLYSIS TO NPC, REMOVE ANY IPSC DATA COLUMNS 
batches=batches[batches$Cell=='NPC',]
data=subset(data,select=c('GeneID',batches$TechRep))

#remove genes with 0 counts 
data=data[rowSums(data[,2:dim(data)[2]])>0,]

data[,2:25]=round(data[,2:25])

#make DESEQ design matrix 
batches$Condition=factor(batches$Condition)
batches$Batch=factor(batches$Batch)

dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData=batches, 
                              design=~Condition+Batch, tidy = TRUE)
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
mat=assay(vsd)
mod1=model.matrix(~Condition,data=batches)

mat=limma::removeBatchEffect(mat,vsd$Batch,design=mod1)
assay(vsd)=mat

dists <- dist(t(assay(vsd)))
plot(hclust(dists))

DESeq2::plotPCA(vsd,intgroup=c("Batch"))
DESeq2::plotPCA(vsd,intgroup=c("Condition"))

data_stabilized=getVarianceStabilizedData(dds)
write.table(mat,file="NPC_only.VST_corrected.txt",sep='\t',row.names = T,col.names = T,quote=F)

saveRDS(vsd, file = "NPC_vsd.rds")