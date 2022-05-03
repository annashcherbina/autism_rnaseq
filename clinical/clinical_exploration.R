rm(list=ls())
#Load latest version of heatmap.3 function
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
library(ggplot2)
library(gplots)
clin_data=read.table("APPiPSC_clinicalMRIdata_031522.csv",header=T,sep='\t')
metrics=c("subj_id",
          "dm_status",
          "visit",
          "ados_severity",
          "ados_sa_severity",
          "ados_rrb_severity",
          "gi_symptoms",
          "gi_severity",
          "vitals_hc_cm",
          "vitals_wt_lb",
          "vitals_bmi",
          "mori_total_l1_cerebellum",
          "mori_total_l1_hemisphere_r",
          "mori_total_l1_hemisphere_l",
          "mori_total_l1_total_volume",
          "mori_total_l1_csf",
          "mori_total_l1_brainstem",
          "mori_total_l3_frontal_l",
          "mori_total_l3_frontal_r",
          "mori_total_l3_parietal_l",
          "mori_total_l3_parietal_r",
          "mori_total_l3_temporal_l",
          "mori_total_l3_temporal_r",
          "mori_total_l3_limbic_l",
          "mori_total_l3_limbic_r",
          "mori_total_l3_occipital_l",
          "mori_total_l3_occipital_r",
          "mori_total_l3_lateralventricle_l",
          "mori_total_l3_lateralventricle_r",
          "mori_total_l3_iii_ventricle",
          "mori_total_l3_iv_ventricle",
          "iq_viq",
          "iq_nviq",
          "iq_fsiq",
          "vine_communication_ss",
          "vine_living_ss",
          "vine_social_ss",
          "vine_motor_ss",
          "vine_abc_ss")
clin_metrics=clin_data[,metrics]
clin_metrics=clin_metrics[order(clin_metrics$dm_status,clin_metrics$subj_id, clin_metrics$visit),]
row.names(clin_metrics)=paste(clin_metrics$subj_id,clin_metrics$visit,sep='-')
#get the disease assignments 
condition=clin_metrics$dm_status

rowsidecolors=data.frame(condition)

#replace rowsidecolors with color names
rowsidecolors[rowsidecolors=="TD-N"]="#FF0000"
rowsidecolors[rowsidecolors=="ASD-N"]="#00FF00"
rowsidecolors[rowsidecolors=="ASD-DM"]="#0000FF"
colnames(rowsidecolors)=c("Condition")

clin_metrics$dm_status=NULL
clin_metrics$subj_id=NULL
clin_metrics$visit=NULL

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")
#normalize metrics across subjects 
require(gtools)
require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

rowsidecolors=as.matrix(rowsidecolors)
clin_metrics=as.matrix(clin_metrics)
heatmap.3(clin_metrics,
          trace="none",
          scale="column",
          Rowv=FALSE,
          Colv=FALSE,
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="",
          RowSideColorsSize = 2,
          RowSideColors = t(as.matrix(rowsidecolors)),
          dendrogram = 'none',
          symbreak=FALSE,
          margins=c(20,10))

legend('left',legend=c("TDN","ASDN","ASDDM"),
       fill=c("#FF0000","#00FF00","#0000FF"),
       border=FALSE, bty="n", y.intersp = 0.8, cex=0.7)

############################################
### FILTERED DATA FOR FIRST, LAST, DELTA ### 
############################################
conditions=read.table('subject_to_condition.txt',header=T,sep='\t')
conditions=conditions[order(conditions$Condition),]
rowsidecolors=data.frame(conditions$Condition)
#replace rowsidecolors with color names
rowsidecolors[rowsidecolors=="TD-N"]="#FF0000"
rowsidecolors[rowsidecolors=="ASD-N"]="#00FF00"
rowsidecolors[rowsidecolors=="ASD-DM"]="#0000FF"
colnames(rowsidecolors)=c("Condition")

first_data=read.table("first_vals.txt",header=T,sep='\t',row.names = 1)
first_data=first_data[conditions$Subject,]

last_data=read.table("last_vals.txt",header=T,sep='\t',row.names=1)
last_data=last_data[conditions$Subject,]

delta_data=read.table("delta_vals.txt",header=T,sep='\t',row.names=1)
common_subjects=intersect(conditions$Subject,
                          row.names(delta_data))
delta_data=delta_data[common_subjects,]
delta_rowsidecolors=as.data.frame(conditions[conditions$Subject %in% common_subjects,c("Condition")])
#replace rowsidecolors with color names
delta_rowsidecolors[delta_rowsidecolors=="TD-N"]="#FF0000"
delta_rowsidecolors[delta_rowsidecolors=="ASD-N"]="#00FF00"
delta_rowsidecolors[delta_rowsidecolors=="ASD-DM"]="#0000FF"
colnames(delta_rowsidecolors)=c("Condition")


heatmap.3(first_data,
          trace="none",
          scale="column",
          Rowv=FALSE,
          Colv=FALSE,
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="First Observation",
          RowSideColorsSize = 2,
          RowSideColors = t(as.matrix(rowsidecolors)),
          dendrogram = 'none',
          symbreak=FALSE,
          margins=c(20,10))

legend('left',legend=c("TDN","ASDN","ASDDM"),
       fill=c("#FF0000","#00FF00","#0000FF"),
       border=FALSE, bty="n", y.intersp = 0.8, cex=0.7)
heatmap.3(last_data,
          trace="none",
          scale="column",
          Rowv=FALSE,
          Colv=FALSE,
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="Final Observation",
          RowSideColorsSize = 2,
          RowSideColors = t(as.matrix(rowsidecolors)),
          dendrogram = 'none',
          symbreak=FALSE,
          margins=c(20,10))

legend('left',legend=c("TDN","ASDN","ASDDM"),
       fill=c("#FF0000","#00FF00","#0000FF"),
       border=FALSE, bty="n", y.intersp = 0.8, cex=0.7)

heatmap.3(delta_data,
          trace="none",
          scale="column",
          Rowv=FALSE,
          Colv=FALSE,
          distfun=distCor,
          hclustfun=hclustAvg,
          col=rev(cols),
          main="Delta Observation",
          RowSideColorsSize = 2,
          RowSideColors = t(as.matrix(delta_rowsidecolors)),
          dendrogram = 'none',
          symbreak=FALSE,
          margins=c(20,10))

legend('left',legend=c("TDN","ASDN","ASDDM"),
       fill=c("#FF0000","#00FF00","#0000FF"),
       border=FALSE, bty="n", y.intersp = 0.8, cex=0.7)
