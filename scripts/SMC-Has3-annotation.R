## Load libraries----
library(Seurat)

## Load Seurat object----
load("results/objects/obj.Rdata")

## Cluster annotation----
DimPlot(obj, label = T)
obj <- RenameIdents(obj,
  `0` = "SMC-Mod1",
  `1` = "SMC-Mod2",
  `2` = "SMC-Mod3",
  `3` = "SMC-Mod4",
  `4` = "SMC-Con",
  `5` = "EC",
  `6` = "MAC")
obj@meta.data$annotation <- Idents(obj)
DimPlot(obj, label = T, group.by = "annotation")
Idents(obj) <- "annotation"

## Save annotated Seurat object----
save(obj, file = "results/objects/obj.Rdata")
