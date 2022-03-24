## Load libraries----
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)

## Load Seurat object----
load("results/objects/obj.Rdata")

## Annotated cluster markers----
idents_list <- levels(obj)
conserved_markers_list <- list()
for (i in idents_list) {
  cluster <- i
  markers <- FindConservedMarkers(obj,
    ident.1 = cluster, min.pct = 0.25, logfc.threshold = 0.25,
    only.pos = T, assay = "RNA", slot = "data", grouping.var = "genotype"
  )
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
  file = "results/cluster-markers/conserved_markers.csv",
  row.names = F
)

## Dotplot of conserved marker genes----

# select bottom combined p value
top10 <- conserved_markers %>%
  group_by(cluster) %>%
  top_n(n = -10, wt = minimump_p_val)

p <- DotPlot(
  object = obj, features = rev(top10$gene),
  cols = c("#4878CD", "#D75438"),
  dot.scale = 6,
  split.by = "genotype"
) +
  RotatedAxis() +
  xlab("Marker genes") +
  theme(
    axis.text.x = element_text(face = "italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  geom_vline(xintercept = c(10.5, 20.5, 30.5, 40.5, 50.5, 60.5))

pdf(
  file = "results/cluster-markers/DotPlot_top10_genotype.pdf",
  height = 5,
  width = 18,
  useDingbats = F
)
print(p)
dev.off()
