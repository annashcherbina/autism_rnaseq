library(EnhancedVolcano)
volcano_wrapper=function(df,genes_to_highlight=c(),title="",
                         pval_thresh=0.01,
                         lfc_thresh=1,
                         xlim=c(-10,10),
                         ylim=c(0,10)){
  library(ggrepel)
  options(ggrepel.max.overlaps=Inf)
  keyvals <- rep('#CCCCCC', nrow(df))
  names(keyvals) <- rep('None', nrow(df))
  keyvals[which((df$logFC > lfc_thresh) &
                  (df$adj.P.Val < pval_thresh))] <- '#CA0020'
  names(keyvals)[which((df$logFC > lfc_thresh) &
                         (df$adj.P.Val < pval_thresh))] <-
    rep('Up', length(which(
      (df$logFC > lfc_thresh) & (df$adj.P.Val < pval_thresh)
    )))
  keyvals[which((df$logFC < -lfc_thresh) &
                  (df$adj.P.Val < pval_thresh))] <- '#0571B0'
  names(keyvals)[which((df$logFC < -lfc_thresh) &
                         (df$adj.P.Val < pval_thresh))] <-
    rep('Down', length(which(
      (df$logFC < -lfc_thresh) & (df$adj.P.Val < pval_thresh)
    )))
  
  EnhancedVolcano(
    df,
    lab = df$Gene,
    drawConnectors=TRUE,
    selectLab=genes_to_highlight,
    x = 'logFC',
    y = 'adj.P.Val',
    pCutoff = pval_thresh,
    FCcutoff = lfc_thresh,
    title = title,
    labCol = 'black',
    labFace = 'bold',
    caption = "",
    colCustom = keyvals,
    
    xlim=xlim,
    ylim=ylim,
    labSize = 4)
}