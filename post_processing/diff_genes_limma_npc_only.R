rm(list=ls())
library(limma)
library(sva)
library(ggplot2)
library(matrixStats)
data=read.table("merged_rsem/rna.counts.txt.tpm",header=TRUE,check.names = F)
batches=read.table("merged_rsem/batches.txt",header=TRUE)

data=aggregate(data[,3:ncol(data)],by=list(GeneID=data$GeneName),FUN=sum)
row.names(data)=data$GeneID
data$GeneID=NULL

batches$TechRep=as.character(batches$TechRep)
data=data[,batches$TechRep[batches$Cell=="NPC"]]
batches=batches[batches$Cell=="NPC",]
row.names(batches)=batches$TechRep

print(dim(data))

#remove genes with TPM < 1 in all samples
data=data[rowMaxs(as.matrix(data))>1,]
print(dim(data))

#remove lncRNA, snoRNA, and similar
genes_to_keep=read.table("gencode_to_category_filtered.txt",header=T)$g_name
data=data[row.names(data) %in% genes_to_keep,]
print(dim(data))


#remove low variance genes
data <- genefilter::varFilter(as.matrix(data), var.func = IQR,
                              var.cutoff = 0.01, filterByQuantile = FALSE)
print(dim(data))

#log-transform
data=log2(data+0.1)


#for cor, replace 0 with NA to avoid correlating 0
data_forcor=data
#data_forcor[data_forcor==0]=NA

pairwise_cor_data_spearman=round(cor(data_forcor,method='spearman',use='pairwise.complete.obs'),2)
pairwise_cor_data_pearson=round(cor(data_forcor,method='pearson',use='pairwise.complete.obs'),2)

library(pheatmap)
batches$Sample=factor(batches$Sample)
pheatmap(pairwise_cor_data_spearman, 
         display_numbers = T,
         cluster_rows = T,
         cluster_cols = T,
         annotation_row=batches[,c("Batch","Condition","Sample")],
         main='Spearman')
pheatmap(pairwise_cor_data_pearson, 
         display_numbers = T,
         cluster_rows = T,
         cluster_cols = T,
         annotation_row=batches[,c("Batch","Condition","Sample")],
         main='Pearson')
data.pca=prcomp(t(data))
var_explained=as.character(round(100*data.pca$sdev^2/sum(data.pca$sdev^2),2))
barplot(100*data.pca$sdev^2/sum(data.pca$sdev^2),las=2,ylab="% Variance Explained",xlab="Principal Component")

pca_df=data.frame(data.pca$x)
pca_df=cbind(pca_df,batches)
pca_df$Sample=factor(pca_df$Sample)
p1=ggplot(data=pca_df,aes(x=pca_df$PC1,
                       y=pca_df$PC2,
                      color=pca_df$Sample))+
  geom_point(size=3)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: PC1 vs PC2")+
  theme_bw()

p2=ggplot(data=pca_df,aes(x=pca_df$PC1,
                          y=pca_df$PC2,
                          color=pca_df$Batch))+
  geom_point(size=3)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: PC1 vs PC2")+
  theme_bw()

p3=ggplot(data=pca_df,aes(x=pca_df$PC1,
                          y=pca_df$PC2,
                          color=pca_df$Condition))+
  geom_point(size=4)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: PC1 vs PC2")+
  theme_bw(20)+scale_color_manual(name='Group',values=c('#377eb8','#e41a1c','#4daf4a'))+
  theme(legend.position = "bottom")

source('~/helpers.R')
multiplot(p1,p2,p3,cols=3)

#DESIGN MATRIX
mod0=model.matrix(~1,data=batches)
mod1=model.matrix(~0+Condition,data=batches)

#RUN SVA 
sva.obj=sva(data,mod1,mod0)
sur_var=data.frame(sva.obj$sv)
cleaned_E=removeBatchEffect(data,covariates=sur_var, batch=batches$Batch,design=mod1)

data.pca=prcomp(t(cleaned_E),center=T,scale=T)
var_explained=as.character(round(100*data.pca$sdev^2/sum(data.pca$sdev^2),2))
barplot(100*data.pca$sdev^2/sum(data.pca$sdev^2),las=2,ylab="% Variance Explained",xlab="Principal Component")
pca_df=data.frame(data.pca$x)
pca_df=cbind(pca_df,batches)
library()
p1=ggplot(data=pca_df,aes(x=pca_df$PC1,
                          y=pca_df$PC2,
                          color=pca_df$Sample))+
  geom_point(size=3)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: PC1 vs PC2")+
  theme_bw()

p2=ggplot(data=pca_df,aes(x=pca_df$PC1,
                          y=pca_df$PC2,
                          color=pca_df$Batch))+
  geom_point(size=3)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA: PC1 vs PC2")+
  theme_bw()

pca_df$ Condition=as.character(pca_df$Condition)
pca_df$Condition[pca_df$Condition=='TDN']='TD-N'
pca_df$Condition[pca_df$Condition=='ASDDM']='ASD-DM'
pca_df$Condition[pca_df$Condition=='ASDN']='ASD-N'


p3=ggplot(data=pca_df,aes(x=pca_df$PC1,
                          y=pca_df$PC2,
                          color=pca_df$Condition))+
  geom_point(size=3)+
  xlab(paste("PC1:",var_explained[1]))+
  ylab(paste("PC2:",var_explained[2]))+
  ggtitle("PCA post SVA correction")+
  scale_color_manual(name="Group",values=c('#377eb8','#4daf4a','#e41a1c'))+
  theme_bw()+
  theme(legend.position="bottom")
ggsave('figure1_L.pdf',p3,width = 4, height = 4, units = "in", device = cairo_pdf())
multiplot(p1,p2,p3,cols=3)
batches$Condition=paste0("NPC",batches$Condition)
#GET VARIANCE ACROSS SAMPLES RELATIVE TO THE MEAN 
tdn_sub=2^cleaned_E[,batches$TechRep[batches$Condition=="NPCTDN"]]
asdn_sub=2^cleaned_E[,batches$TechRep[batches$Condition=="NPCASDN"]]
asddm_sub=2^cleaned_E[,batches$TechRep[batches$Condition=="NPCASDDM"]]

tdn_sub_var_mean_ratio=rowVars(tdn_sub)/rowMeans(tdn_sub)
asdn_sub_var_mean_ratio=rowVars(asdn_sub)/rowMeans(asdn_sub)
asddm_sub_var_mean_ratio=rowVars(asddm_sub)/rowMeans(asddm_sub)


#LIMMA REGRESSION ANALYSIS
fit=lmFit(cleaned_E,mod1)

source('volcano_wrapper.R')

pthresh=0.01
lfc_thresh=1 #fold change of 2
mod1=model.matrix(~0+Condition,data=batches)

g1=c("ConditionNPCASDDM","ConditionNPCASDDM","ConditionNPCASDN")
g2=c("ConditionNPCASDN","ConditionNPCTDN","ConditionNPCTDN")
for(i in seq(1,3))
{
  comparison=paste(g1[i],'vs',g2[i],sep='_')
  print(comparison)
  v=rep(0,ncol(mod1))
  v[which(colnames(mod1)==g1[i])] = 1
  v[which(colnames(mod1)==g2[i])]=-1
  contr = as.matrix(v)
  colnames(contr) = "Group"
  rownames(contr) = colnames(mod1)
  fit0 <- limma::contrasts.fit(fit, contrasts = contr)
  fit0 <- limma::eBayes(fit0)#, trend = TRUE, proportion = 0.05)
  de = limma::topTable(
    fit = fit0,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none",
  )
  de$TDN_var_mean_ratio=tdn_sub_var_mean_ratio
  de$ASDN_var_mean_ratio=asdn_sub_var_mean_ratio
  de$ASDDM_var_mean_ratio=asddm_sub_var_mean_ratio
  
  de=na.omit(de)
  de$abst=abs(de$t)
  de$Gene=row.names(de)
  #get the DE sig genes
  de_sig=de[(de$adj.P.Val<pthresh) & (abs(de$logFC)>lfc_thresh), ]
  #write to output file
  write.table(de_sig,file=paste(comparison,'.de.sig.tsv',sep=''),sep='\t',row.names=T,col.names=T)
  write.table(de,file=paste(comparison,'.de.tsv',sep=''),sep='\t',row.names=T,col.names=T)
  
  #make volcano plot
  de_sig=de_sig[order(de_sig$abst,decreasing = T),]
  labels=row.names(de_sig)[1:20]
  up_count=paste('up',sum(as.integer(de_sig$logFC>0)))
  down_count=paste('down',sum(as.integer(de_sig$logFC<0)))
  title=paste(comparison,up_count,down_count,sep='\n')
  
  print(volcano_wrapper(de,
                        title=title,
                        genes_to_highlight = labels,
                        pval_thresh=pthresh,
                        lfc_thresh=lfc_thresh,
                        xlim=c(-8,8),
                        ylim=c(0,15)))
  
  
  #get pathway enrichments
  foreground=de_sig$Gene
  background=de$Gene
  out=gProfileR::gprofiler(query = foreground,
                           custom_bg = background,
                           exclude_iea = T,
                           src_filter = c("GO:BP"),
                           organism = "hsapiens",
                           hier_filtering = "moderate")
  out$logp=-1*log10(out$p.value)
  out=out[order(out$p.value,decreasing = T),]
  out$term.name=factor(out$term.name,levels=out$term.name)
  out=out[out$term.size<500,]
  write.table(out,paste(comparison,'.gProfiler.tsv',sep=''),sep='\t',col.names=T,row.names=T)
  #plot pathway enrichments
  print(ggplot(out,aes(x=term.name,y=logp))+geom_bar(stat='identity')+coord_flip()+xlab("Enriched GO Term")+ylab("-log10(pvalue")+theme_bw(20))
  vis_go=out[,c('term.size','p.value','term.name','term.id','intersection')]
  print(vis_go)
}

write.table(2^cleaned_E,file="NPC_only.corrected_tpm.txt",row.names=TRUE,col.names=TRUE,sep='\t')