library(ggplot2)
genes_of_interest=c("ENSG00000138646.9-HERC5",
                    "ENSG00000119917.14-IFIT3",
                    "ENSG00000185201.16-IFITM2",
                    "ENSG00000187608.10-ISG15",
                    "ENSG00000172183.15-ISG20",
                    "ENSG00000157601.14-MX1",
                    "ENSG00000002586.20-CD99",
                    "ENSG00000079616.12-KIF22",
                    "ENSG00000103495.14-MAZ",
                    "ENSG00000124181.14-PLCG1",
                    "ENSG00000136068.14-FLNB",
                    "ENSG00000149923.14-PPP4C",
                    "ENSG00000149927.18-DOC2A",
                    "ENSG00000167371.21-PRRT2",
                    "ENSG00000174938.14-SEZ6L2")
genes_of_interest=c("ENSG00000120885.22-CLU",
                    "ENSG00000137331.12-IER3",
                    "ENSG00000067225.18-PKM",
                    "ENSG00000135549.15-PKIB")
data=read.table("merged_rsem/rna.counts.txt.tpm",header=TRUE,sep='\t')
# columns to paste together
cols <- c( 'GeneID' , 'GeneName')

# create a new column `x` with the three columns collapsed together
rownames(data) <- apply( data[ , cols ] , 1 , paste , collapse = "-" )

# remove the unnecessary columns
data <- data[ , !( names( data ) %in% cols ) ]

data=data[genes_of_interest,]
data=as.data.frame(t(data))
batches=read.table("merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1)

data=merge(data,batches,by=0)

conditions = c("TDN","ASDN","ASDDM")

data$Condition=factor(data$Condition, levels=conditions)

celltype = c("IPSC","NPC")

plot.list = sapply(genes_of_interest, function(gene){
  p.list = lapply(celltype, function(ct){
    
    this_data = data[data$Cell==ct, , drop = F]
    p <- ggplot(data=this_data,
                aes(x=this_data$Condition,
                    y=this_data[[gene]]))+
      geom_boxplot(color='#FF0000',outlier.shape=NA)+
      geom_jitter()+
      ylab("TPM")+
      xlab("Condition")+
      ggtitle(paste(gene,'\nTPM,',ct))+
      theme(plot.title = element_text(size = 8))
    return(p)
  })
  return(p.list)
}, simplify = T)

library(gridExtra)
grid.arrange(grobs = plot.list,ncol=2)
