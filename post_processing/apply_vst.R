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
batches$id=batches$TechRep
dds <- DESeqDataSetFromMatrix(countData=data, 
                              colData=batches, 
                              design=~Condition, tidy = TRUE)

vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
dists <- dist(t(assay(vsd)))
plot(hclust(dists))


dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)


data_stabilized=getVarianceStabilizedData(dds)
write.table(data_stabilized,file="NPC_only.VST_corrected.txt",sep='\t',row.names = T,col.names = T,quote=F)
