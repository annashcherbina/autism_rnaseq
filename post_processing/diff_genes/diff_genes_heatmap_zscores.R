library(ggplot2)
library(gplots)
library("gplots")
library("devtools")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

source("~/helpers.R")

topn=40
asdn_vs_tdn=read.table('ConditionASDN_vs_ConditionTDN.de.sig.varfilter.tsv',header=T,sep='\t',check.names = F)
asdn_vs_tdn=asdn_vs_tdn[order(asdn_vs_tdn$adj.P.Val),][1:topn,]
asddm_vs_tdn=read.table('ConditionASDDM_vs_ConditionTDN.de.sig.varfilter.tsv',header=T,sep='\t',check.names = F)
asddm_vs_tdn=asddm_vs_tdn[order(asddm_vs_tdn$adj.P.Val),][1:topn,]

asddm_vs_asdn=read.table('ConditionASDDM_vs_ConditionASDN.de.sig.varfilter.tsv',header=T,sep='\t',check.names=F)
asddm_vs_asdn=asddm_vs_asdn[order(asddm_vs_asdn$adj.P.Val),][1:topn,]

top_sig_genes=unique(c(asddm_vs_tdn$Gene,asddm_vs_asdn$Gene, asdn_vs_tdn$Gene))
data=read.table("../NPC_only.corrected_tpm.txt",header=T,sep='\t',row.names = 1, check.names = F)[top_sig_genes,]

batches=read.table("../merged_rsem/batches.txt",header=TRUE,sep='\t',check.names = F)
batches=batches[batches$Cell=="NPC",]
batches=batches[order(batches$Condition),]

data=data[,batches$TechRep]

colsidecolors=data.frame(batches$Condition)

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
heatmap.3(data,
          trace="none",
          scale="row",
          Rowv=TRUE,
          Colv=FALSE,
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="",
          ColSideColorsSize = 1,
          ColSideColors = as.matrix(colsidecolors),
          symbreak=FALSE,
          cexRow=0.7,
          margins=c(5,20))

legend('topright',legend=c("TDN","ASDN","ASDDM"),
       fill=c("#FF0000","#00FF00","#0000FF"),
       border=FALSE, bty="n", y.intersp = 0.6, cex=0.7)
#dev.off()