## Load libraries----
library(Seurat)
library(ggplot2)

## Load Seurat object----
load("results/objects/obj.Rdata")

## Dimplot----
dimcols <- c("#EDAE3C",
             "#D95030",
             "#FB5938",
             "#D6809D",
             "#508AB2",
             "#6FA7A0",
             "#7572AD")
p <- DimPlot(obj,
                label = F,
                pt.size = 1,
                cols = dimcols) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  NoLegend()
p <- LabelClusters(plot = p,
              id = "ident",
              repel = T,
              force = 10,
              size = 6)
pdf(file = "results/dimplots/dimplot.pdf",
    width = 6,
    height = 6,
    useDingbats = F)
print(p)
dev.off()

## DimPlot grouped by genotype----
obj@meta.data$genotype <- factor(obj@meta.data$genotype, levels = c("WT", "KO"))
p <- DimPlot(obj,
             group.by = "genotype",
             pt.size = 1,
             cols = c("#4878CD", "#D75438")) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 20),
        legend.position = c(0.075, 0.15),
        legend.text = element_text(size = 20),
        plot.title = element_blank())
pdf(file = "results/dimplots/dimplot_genotype.pdf",
    width = 6,
    height = 6,
    pointsize = 12,
    useDingbats = F)
print(p)
dev.off()

## Save factored Seurat object----
save(obj, file = "results/objects/obj.Rdata")