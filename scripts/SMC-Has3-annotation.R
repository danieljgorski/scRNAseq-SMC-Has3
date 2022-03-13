## Load libraries----
library(Seurat)
library(multtest)

## Load Seurat object----
load("results/objects/obj.Rdata")

## List of all genes----
all_genes <- rownames(obj)
write.table(all_genes,
            file = "results/cluster_markers/all_genes.txt",
            row.names = F,
            col.names = F,
            quote = F)

## Find conserved cluster markers----
Idents(obj) <- "seurat_clusters"
DefaultAssay(obj) <- "RNA"
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
conserved_markers <- do.call(rbind, conserved_markers_list)
remove(conserved_markers_list)
write.csv(conserved_markers, file = "results/cluster_markers/conserved_markers.csv", row.names = F)

## Cluster annotation----
DimPlot(obj, label = T)
obj <- RenameIdents(obj,
                    "0" = "SMC-Mod1",
                    "1" = "SMC-Mod2",
                    "2" = "SMC-Mod3",
                    "3" = "SMC-Mod4",
                    "4" = "SMC-Con",
                    "5" = "EC",
                    "6" = "MAC")
obj@meta.data$annotation <- Idents(obj)
DimPlot(obj, label = T, group.by = "annotation")
Idents(obj) <- "annotation"

## Save annotated Seurat object----
save(obj, file = "results/objects/obj.Rdata")