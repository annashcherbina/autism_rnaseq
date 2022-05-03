rm(list=ls())
library(ggplot2)
library(gProfileR)
data=read.table("../post_processing/merged_rsem/rna.counts.txt.tpm",header=TRUE,check.names = F)
batches=read.table("../post_processing/merged_rsem/batches.txt",header=TRUE)
batches$Sample=factor(batches$Sample)
#sum tpm across genes w/ same name
data=aggregate(data[,3:ncol(data)],by=list(GeneID=data$GeneName),FUN=sum)
row.names(data)=data$GeneID
data$GeneID=NULL

#WE RESTRICT OUT ANALLYSIS TO NPC, REMOVE ANY IPSC DATA COLUMNS 
batches=batches[batches$Cell=='NPC',]
batches=batches[order(batches$Condition),]
data=subset(data,select=c(batches$TechRep))

#remove genes with 0 counts 
data=data[rowSums(data)>0,]

#remove lncRNA, snoRNA, and similar
genes_to_keep=read.table("../gencode_to_category_filtered.txt",header=T)$g_name
data=data[row.names(data) %in% genes_to_keep,]

sig_ASDDM_ASDN=read.table('all.sig.ASDDM_ASDN.tsv',header=T,sep='\t')
genes_ASDDM_ASDN=unique(sig_ASDDM_ASDN$Gene)
ASDDM_ASDN=data[genes_ASDDM_ASDN,]

sig_ASDN_TDN=read.table('all.sig.ASDN_TDN.tsv',header=T,sep='\t') 
genes_ASDN_TDN=unique(sig_ASDN_TDN$Gene)
ASDN_TDN=data[genes_ASDN_TDN,]

sig_ASDDM_TDN=read.table('all.sig.ASDDM_TDN.tsv',header=T,sep='\t')
genes_ASDDM_TDN=unique(sig_ASDDM_TDN$Gene)
ASDDM_TDN=data[genes_ASDDM_TDN,]


library("gplots")
library("devtools")
require(gtools)
require(RColorBrewer)

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
colsidecolors=data.frame(batches$Condition)

#replace colsidecolors with color names
colsidecolors[colsidecolors=="TDN"]="#FF0000"
colsidecolors[colsidecolors=="ASDN"]="#00FF00"
colsidecolors[colsidecolors=="ASDDM"]="#0000FF"

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")
heatmap.3(ASDDM_TDN,
          trace="none",
          scale="row",
          dendro='none',
          Rowv=TRUE,
          Colv=FALSE,
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="ASDDM vs TDN",
          ColSideColorsSize = 2,
          ColSideColors = as.matrix(colsidecolors),
          symbreak=FALSE,
          keysize=0,
          margins=c(5,20))

legend('topright',legend=c("TDN","ASDN","ASDDM"),
       fill=c("#FF0000","#00FF00","#0000FF"),
       border=FALSE, bty="n", y.intersp = 0.6, cex=0.7)


heatmap.3(ASDDM_ASDN,
          trace="none",
          scale="row",
          dendro='none',
          Rowv=TRUE,
          Colv=FALSE,
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="ASDDM vs ASDN",
          ColSideColorsSize = 2,
          keysize=0,
          ColSideColors = as.matrix(colsidecolors),
          symbreak=FALSE,
          margins=c(5,20))

legend('topright',legend=c("TDN","ASDN","ASDDM"),
       fill=c("#FF0000","#00FF00","#0000FF"),
       border=FALSE, bty="n", y.intersp = 0.6, cex=0.7)


heatmap.3(ASDN_TDN,
          trace="none",
          scale="row",
          dendro='none',
          Rowv=TRUE,
          Colv=FALSE,
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="ASDN vs TDN",
          ColSideColors = as.matrix(colsidecolors),
          symbreak=FALSE,
          margins=c(5,20))

legend('topright',legend=c("TDN","ASDN","ASDDM"),
       fill=c("#FF0000","#00FF00","#0000FF"),
       border=FALSE, bty="n", y.intersp = 0.6, cex=0.7)


foreground=genes_ASDDM_ASDN
background=row.names(data)

out=gProfileR::gprofiler(query = foreground,
                         custom_bg = background,
                         exclude_iea = T,
                         src_filter = c("GO:BP"),
                         organism = "hsapiens",
                         hier_filtering = "moderate")
out=out[order(out$p.value,decreasing = T),]
out=out[out$term.size<500,]
out$query.number=NULL
out$significant=NULL
out$subgraph.number=NULL



foreground=genes_ASDDM_TDN
background=row.names(data)

out=gProfileR::gprofiler(query = foreground,
                         custom_bg = background,
                         exclude_iea = T,
                         src_filter = c("GO:BP"),
                         organism = "hsapiens",
                         hier_filtering = "moderate")
out=out[order(out$p.value,decreasing = T),]
out=out[out$term.size<500,]
out$query.number=NULL
out$significant=NULL
out$subgraph.number=NULL




foreground=genes_ASDN_TDN
background=row.names(data)

out=gProfileR::gprofiler(query = foreground,
                         custom_bg = background,
                         exclude_iea = T,
                         src_filter = c("GO:BP"),
                         organism = "hsapiens",
                         hier_filtering = "moderate")
out=out[order(out$p.value,decreasing = T),]
out=out[out$term.size<500,]
out$query.number=NULL
out$significant=NULL
out$subgraph.number=NULL