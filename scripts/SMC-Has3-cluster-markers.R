## Load libraries----
library(Seurat)
library(ggplot2)

## Load Seurat object----
load("results/objects/obj.Rdata")

## Annotated cluster markers----
idents_list <- levels(obj)
conserved_markers_list <- list()
for (i in idents_list) {
  cluster <- i
  conserved_markers <- FindConservedMarkers(obj,
                                            ident.1 = cluster,
                                            min.pct = 0.25,
                                            logfc.threshold = 0.25,
                                            only.pos = T,
                                            assay = "RNA",
                                            slot = "data",
                                            grouping.var = "genotype")
  conserved_markers <- conserved_markers[conserved_markers$WT_p_val_adj < 0.05, ]
  conserved_markers <- conserved_markers[conserved_markers$KO_p_val_adj < 0.05, ]
  conserved_markers$cluster <- i
  conserved_markers$gene <- row.names(conserved_markers)
  row.names(conserved_markers) <- NULL
  conserved_markers_list[[i]] <- conserved_markers
  remove(conserved_markers)
}
conserved_markers_annotated <- do.call(rbind, conserved_markers_list)
remove(conserved_markers_list)
write.csv(conserved_markers_annotated,
          file = "results/cluster_markers/conserved_markers_annotated.csv",
          row.names = F)