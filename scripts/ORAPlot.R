ORAPlot <- function(genes_of_interest, background_genes, plot_title, num_categories){
  ora <- enrichGO(gene = genes_of_interest,
                  universe = background_genes,
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  keyType = "SYMBOL",
                  pAdjustMethod = "bonferroni",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = F)
  p <- dotplot(ora,
               showCategory = num_categories,
               label_format = 30,
               orderBy = "x") +
    labs(title = plot_title) +
    theme(legend.position = "right",
          legend.box = "vertical")
  print(p)
}
  