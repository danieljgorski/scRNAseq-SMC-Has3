## Load libraries and functions----
library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)
source("scripts/VolcanoPlot.R")

## Load Seurat object----
load("results/objects/obj.Rdata")

## Volcano plots----

# all modulated SMC KO v WT
dge_mod_genotype <- read_csv("results/differential-gene-expression/dge_mod_genotype.csv")
p <- VolcanoPlot(dge_mod_genotype, "All Modulated SMC: KO vs WT")
pdf(file = "results/volcano-plots/VolcanoPlot_mod_genotype.pdf",
    height = 6.5,
    width = 6.5,
    pointsize = 12,
    useDingbats = F)
print(p)
dev.off()

# modulated v contractile 
dge_mod_con <- read_csv("results/differential-gene-expression/dge_mod_con.csv")
p <- VolcanoPlot(dge_mod_con, "Modulated SMC vs Contractile SMC")
pdf(file = "results/volcano-plots/VolcanoPlot_mod_con.pdf",
    height = 6.5,
    width = 6.5,
    pointsize = 12,
    useDingbats = F)
print(p)
dev.off()

# each cluster KO v WT
dge_genotype <- read_csv("results/differential-gene-expression/dge_genotype.csv")
for (i in levels(Idents(obj))) {
  dge_genotype_cluster <- dge_genotype[dge_genotype$cluster == i,]
  title_cluster <- paste0(i, ": KO vs WT")
  pdf(file = paste0("results/volcano-plots/VolcanoPlot_", i, ".pdf"),
      height = 6.5,
      width = 6.5,
      pointsize = 12,
      useDingbats = F)
  VolcanoPlot(dge_genotype_cluster, title_cluster)
  dev.off()
}
