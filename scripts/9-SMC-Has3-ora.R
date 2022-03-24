## Load libraries and functions----
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(readr)
source("scripts/ORAPlot.R")

## Load Seurat object----
load("results/objects/obj.Rdata")

## GOBP over-representation analysis, cluster markers----
conserved_markers <- read.csv("results/cluster-markers/conserved_markers.csv")
for (i in levels(obj)) {
  goi <- conserved_markers[conserved_markers$cluster == i, ]$gene
  plot_title_cluster <- paste0(i, ": GOBP-ORA")
  pdf(file = paste0("results/ora/ORAPlot_", i, "_markers.pdf"),
    width = 6,
    height = 8,
    useDingbats = F)
  ORAPlot(goi, rownames(obj), plot_title_cluster, 20)
  dev.off()
  }

## GOBP over-representation analysis, modulated SMC KO v WT----

# import dge and filter for significant upregulated genes
dge_mod_genotype <- read_csv("results/differential-gene-expression/dge_mod_genotype.csv")
goi <- dge_mod_genotype[dge_mod_genotype$regulation == "Up",]
goi <- goi[goi$p_val_adj < 0.05,]$gene

# plot
p <- ORAPlot(goi,
              rownames(obj),
              "Upregulated Genes in KO Modulated SMC: GOBP-ORA",
              50)
pdf(file = paste0("results/ora/ORAPlot_mod_genotype.pdf"),
    width = 8,
    height = 10,
    useDingbats = F)
print(p)
dev.off()

## GOBP over-representation analysis, modulated v contractile SMC----

# import dge and filter for significant upregulated genes
dge_mod_con <- read_csv("results/differential-gene-expression/dge_mod_con.csv")
goi <- dge_mod_con[dge_mod_con$regulation == "Up",]
goi <- goi[goi$p_val_adj < 0.05,]$gene

# plot
p <- ORAPlot(goi,
              rownames(obj),
              "Upregulated Genes in Modulated SMC: GOBP-ORA",
              100)
pdf(file = paste0("results/ora/ORAPlot_mod_con.pdf"),
    width = 8,
    height = 15,
    useDingbats = F)
print(p)
dev.off()

