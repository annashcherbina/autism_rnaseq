library(EnhancedVolcano)
volcano_wrapper=function(df,title="",pval_thresh=0.01,lfc_thresh=0.1){
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
  x = 'logFC',
  y = 'adj.P.Val',
  pCutoff = pval_thresh,
  FCcutoff = lfc_thresh,
  xlim = c(-3, 3),
  title = title,
  labCol = 'black',
  labFace = 'bold',
  caption = "",
  colCustom = keyvals,
  labSize = 2
)
}