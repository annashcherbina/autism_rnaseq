library(ggplot2)
library(gplots)
library("gplots")
library("devtools")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

source("~/helpers.R")

data=read.table("asddm_asdn_unique_in_pink_magenta_tan.tpm.txt",header=TRUE,sep='\t',check.names = F)
row.names(data)=data$Gene
data$ENSGID=NULL
data$Gene=NULL

batches=read.table("../merged_rsem/batches.txt",header=TRUE,sep='\t',row.names=1,check.names = F)
batches=batches[batches$Cell=="NPC",]
merged=merge(t(data),batches,by=0)
merged$Condition=factor(merged$Condition,levels=c("TDN","ASDN","ASDDM"))

# Sort by vector name [z] then [x]
merged=merged[
  with(merged, order(Condition, Sample)),
]
condition=as.character(merged$Condition)
merged$Cell=NULL 
merged$Condition=NULL
merged$Batch=NULL
merged$Sample=NULL
merged=as.data.frame(merged)
row.names(merged)=merged$Row.names
merged$Row.names=NULL 
merged=as.matrix(t(merged))
colsidecolors=data.frame(condition)

#replace colsidecolors with color names 
colsidecolors[colsidecolors=="TDN"]="#FF0000"
colsidecolors[colsidecolors=="ASDN"]="#00FF00"
colsidecolors[colsidecolors=="ASDDM"]="#0000FF"
colnames(colsidecolors)=c("Condition")

require(gtools)
require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")
#png(file="genes.NPC.pval.lt.5e-11.fc.7.png",width=10,height=9,units='in',res=300)
png(file="genes.NPC.ASDDM_vs_ASDN.unique.pink.magenta.tan.png",width=5,height=25,units='in',res=160)
heatmap.3(merged,
          trace="none",
          scale="row",
          Rowv=TRUE,
          Colv=FALSE,
          dendrogram = "none",
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="",
          ColSideColorsSize = 2,
          ColSideColors = as.matrix(colsidecolors),
          symbreak=FALSE,
          cexRow=0.8,
          cexCol= 0.8,
          keysize = 0.1,
          margins=c(5,5))
dev.off()






