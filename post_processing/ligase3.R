data=read.table("merged_rsem/tpm.txt",sep='\t',header=T)
data=data[data$Symbol %in% c("ENSG00000005156.12"),]
row.names(data)=data$transcript_id
data$gene_id=NULL
data$transcript_id=NULL
data$Symbol=NULL
data=colSums(data)
