library(gProfileR)
data=read.table('module.enrichments.clinical.assoc.genes.sig.tsv',header=T,sep='\t')
out_vals=list()
print(nrow(data))
for(i in seq(1,nrow(data)))
{
  print(i)
  cur_genes=as.vector(strsplit(data$CommonGenes[i],','))[[1]]
  
  out=gProfileR::gprofiler(query = cur_genes,
                           exclude_iea = T,
                           src_filter = c("GO:BP"),
                           organism = "hsapiens",
                           hier_filtering = "moderate")
  out=out[out$term.size<500,]
  if(nrow(out)>=1){
  out$Comparison=data$Comparison[i] 
  out$ClinFeat=data$ClinFeat[i]
  out$Module=data$Module[i]
  out$Module.P.val=data$P.val[i]
  out$CommonGenes=data$CommonGenes[i]
  out_vals[[i]]=out
  }
}
all_out_vals=do.call('rbind',out_vals)
write.table(all_out_vals,file='gProfiler.for.clinical.correlated.genes.with.module.enrichments.tsv',row.names=T,col.names=T,sep='\t')
