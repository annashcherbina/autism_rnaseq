library(ggplot2)
library(gplots)
library("gplots")
library("devtools")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

source("~/helpers.R")
data=read.table("corrected.tpm.pval.5e-11.lfc7.npc.txt.genes.txt",header=TRUE,sep='\t',row.names = 1,check.names=F)
row.names(data)=data$Gene
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
png(file="genes.NPC.pval.lt.5e-11.fc.7.png",width=10,height=9,units='in',res=300)
heatmap.3(merged,
          trace="none",
          scale="row",
          Rowv=TRUE,
          Colv=FALSE,
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="",
          ColSideColorsSize = 2,
          ColSideColors = as.matrix(colsidecolors),
          symbreak=FALSE,
          cexRow=1,
          margins=c(5,20))

legend('topright',legend=c("TDN","ASDN","ASDDM"),
       fill=c("#FF0000","#00FF00","#0000FF"),
       border=FALSE, bty="n", y.intersp = 0.6, cex=0.7)
dev.off()






