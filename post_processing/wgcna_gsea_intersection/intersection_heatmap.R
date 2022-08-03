library(ComplexHeatmap)
library(ggplot2)
library(gplots)
library("gplots")
library("devtools")
data=read.table("module.gsea.enrichments.filtered.tsv",header=T,sep='\t')
#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

source("~/helpers.R")
data$P.val=NULL

#make a heatmap for each comparison 
asdn_tdn=data[data$Comparison=="ASDN_TDN",]

asddm_asdn=data[data$Comparison=="ASDDM_ASDN",]

asddm_tdn=data[data$Comparison=="ASDDM_TDN",]


asdn_tdn_cast=dcast(asdn_tdn,GOTerm~Module,value.var="Bonferroni",sum)
row.names(asdn_tdn_cast)=asdn_tdn_cast$GOTerm
asdn_tdn_cast$GOTerm=NULL

asddm_asdn_cast=dcast(asddm_asdn,GOTerm~Module,value.var="Bonferroni",sum)
row.names(asddm_asdn_cast)=asddm_asdn_cast$GOTerm
asddm_asdn_cast$GOTerm=NULL


asddm_tdn_cast=dcast(asddm_tdn,GOTerm~Module,value.var="Bonferroni",sum)
row.names(asddm_tdn_cast)=asddm_tdn_cast$GOTerm
asddm_tdn_cast$GOTerm=NULL


require(gtools)
require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

heatmap.3(as.matrix(asdn_tdn_cast),
          trace="none",
          scale="none",
          Rowv=FALSE,
          Colv=FALSE,
          col=rev(cols),
          main="ASDN_TDN",
          symbreak=FALSE,
          cexRow = 0.8,
          cexCol = 1,
          margins=c(5,50),
          )
legend('topright',border=FALSE, bty="n", y.intersp = 0.6, cex=0.7)


library(pheatmap)   
library(gplots)
if (nrow(table) > 100) stop("Too many rows for heatmap, who can read?!")
table=as.matrix(asdn_tdn_cast)
fontsize_row = 10 - nrow(table) / 15
pheatmap(table, col=bluered(256), main="ASDN_TDN", cluster_cols=F, cluster_rows = F,
         fontsize_row=fontsize_row, border_color=NA)

table=as.matrix(asddm_asdn_cast)
max_vals= apply(table, 1, FUN = max)
table=table[max_vals>10,]
fontsize_row = 10 - nrow(table) / 15
pheatmap(table, col=bluered(256), main="ASDDM_ASDN", cluster_cols=F, cluster_rows = F,
         fontsize_row=fontsize_row, border_color=NA)

table=as.matrix(asddm_tdn_cast)
max_vals= apply(table, 1, FUN = max)
table=table[max_vals>5,]

fontsize_row = 10 - nrow(table) / 15
pheatmap(table, col=bluered(256), main="ASDDM_TDN", cluster_cols=F, cluster_rows = F,
         fontsize_row=fontsize_row, border_color=NA)

turquoise=data[data$Module=="turquoise",]
turquoise_cast=dcast(turquoise,GOTerm~Comparison,value.var="Bonferroni",sum)
turquoise_cast=turquoise_cast[order(turquoise_cast$GOTerm),]
row.names(turquoise_cast)=turquoise_cast$GOTerm
turquoise_cast$GOTerm=NULL
table=as.matrix(turquoise_cast)
max_vals= apply(table, 1, FUN = max)
table=table[max_vals>5,]
fontsize_row = 10 - nrow(table) / 15
pheatmap(table, col=bluered(256), main="Turquoise GO Terms", cluster_cols=F, cluster_rows = F,
         fontsize_row=fontsize_row, border_color=NA)


pink=data[data$Module=="pink",]
pink_cast=dcast(pink,GOTerm~Comparison,value.var="Bonferroni",sum)
pink_cast=pink_cast[order(pink_cast$GOTerm),]
row.names(pink_cast)=pink_cast$GOTerm
pink_cast$GOTerm=NULL
table=as.matrix(pink_cast)
max_vals= apply(table, 1, FUN = max)
fontsize_row = 10 - nrow(table) / 15
pheatmap(table, col=bluered(256), main="Pink GO Terms", cluster_cols=F, cluster_rows = F,
         fontsize_row=fontsize_row, border_color=NA)



yellow=data[data$Module=="yellow",]
yellow_cast=dcast(yellow,GOTerm~Comparison,value.var="Bonferroni",sum)
yellow_cast=yellow_cast[order(yellow_cast$GOTerm),]
row.names(yellow_cast)=yellow_cast$GOTerm
yellow_cast$GOTerm=NULL
table=as.matrix(yellow_cast)
max_vals= apply(table, 1, FUN = max)
fontsize_row = 10 - nrow(table) / 15
pheatmap(table, col=bluered(256), main="Yellow GO Terms", cluster_cols=F, cluster_rows = F,
         fontsize_row=fontsize_row, border_color=NA)

