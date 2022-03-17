## Load libraries----
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

## Load Seurat object----
load("results/objects/obj.Rdata")

## Annotated cluster markers----
idents_list <- levels(obj)
conserved_markers_list <- list()
for (i in idents_list) {
  cluster <- i
  markers <- FindConservedMarkers(obj,
                                  ident.1 = cluster,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25,
                                  only.pos = T,
                                  assay = "RNA",
                                  slot = "data",
                                  grouping.var = "genotype")
  markers <- markers[markers$WT_p_val_adj < 0.05, ]
  markers <- markers[markers$KO_p_val_adj < 0.05, ]
  markers$cluster <- i
  markers$gene <- row.names(markers)
  row.names(markers) <- NULL
  conserved_markers_list[[i]] <- markers
  remove(markers)
}
conserved_markers <- do.call(rbind, conserved_markers_list)
remove(conserved_markers_list)
write.csv(conserved_markers,
          file = "results/cluster_markers/conserved_markers.csv",
          row.names = F)

## Gene ontology BP over-representation analysis of cluster markers----
for (i in levels(obj)) {
  goi <- conserved_markers[conserved_markers$cluster == i, ]$gene
  ora <- enrichGO(gene =  goi,
                  universe = rownames(obj),
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  keyType = "SYMBOL",
                  pAdjustMethod = "bonferroni",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable = F)
  p <- dotplot(ora,
          showCategory = 20,
          label_format = 30,
          orderBy = "x") +
  labs(title = i, " marker gene ORA") +
  theme(legend.position = "right",
        legend.box = "vertical")
  pdf(file = paste0("results/cluster_markers/", i, "_marker_gene_ORA.pdf"),
      width = 6,
      height = 8,
      useDingbats = F)
  print(p)
  dev.off()
}