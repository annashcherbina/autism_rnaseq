rm(list=ls())
library(rrvgo)
library(ggplot2)

scatterPlot=function(simMatrix, reducedTerms, size = "score", addLabel = TRUE, 
                     labelSize = 3) 
{
  x <- cmdscale(as.matrix(as.dist(1 - simMatrix)), eig = TRUE, 
                k = 2)
  df <- cbind(as.data.frame(x$points), reducedTerms[match(rownames(x$points), 
                                                          reducedTerms$go), c("term", "parent", "parentTerm", "size")])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = V1, y = V2, color = parentTerm)) + 
    ggplot2::geom_point(ggplot2::aes(size = size), alpha = 0.5) + 
    ggplot2::scale_color_discrete(guide='none') + ggplot2::scale_size_continuous(guide = "none", 
                                                                                 range = c(0, 25)) + ggplot2::scale_x_continuous(name = "") + 
    ggplot2::scale_y_continuous(name = "") + ggplot2::theme_minimal() + 
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                   axis.text.y = ggplot2::element_blank())
  if (addLabel) {
    p + ggrepel::geom_label_repel(aes(label = parentTerm), 
                                  max.overlaps = 20,
                                  data = subset(df, parent == rownames(df)), box.padding = grid::unit(1, 
                                                                                                      "lines"), size = labelSize)
  }
  else {
    p
  }
}



samples=c("ConditionASDDM_vs_ConditionASDN","ConditionASDN_vs_ConditionTDN","ConditionASDDM_vs_ConditionTDN")
for(sample in samples)
{
gsea_data=read.table(paste0('revigo/',sample,'.gsea.tsv'),header=T,sep='\t')
  simMatrix <- calculateSimMatrix(gsea_data$GOterm,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

  scores <- setNames(-log10(gsea_data$enrichment_P.value), gsea_data$GOterm)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

  write.table(reducedTerms,paste0(sample,".clustered.go.tsv"),quote = F,sep='\t',row.names = F)


  png(paste0(sample,".clustered.go.scatter.png"))
  scatterPlot(simMatrix, reducedTerms)
  dev.off()
  
  #png(paste0(sample,".clustered.go.treemap.png"))
  #treemapPlot(reducedTerms)
  #dev.off()
}