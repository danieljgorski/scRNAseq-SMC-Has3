## Load libraries----
library(Seurat)
library(ggplot2)
library(dplyr)

## Load Seurat object----
load("results/objects/obj.Rdata")

## DGE KO v WT----

# initiate a list to hold results
dge <- list()

# loop through each cluster, find DGE with no threshold 
for (i in levels(Idents(obj))) {
  results <- FindMarkers(obj,
                         subset.ident = i,
                         group.by = "genotype",
                         ident.1 = "KO",
                         assay = "RNA",
                         slot = "data",
                         logfc.threshold = 0,
                         verbose = T)
  results$cluster <- i
  results$gene <- row.names(results)
  row.names(results) <- NULL
  results$regulation <- ifelse(results$avg_logFC > 0, "Up", "Down")
  dge[[i]] <- results
}

# row bind the list
dge_genotype <- do.call(rbind, dge)
rownames(dge_genotype) <- NULL

# save
write.csv(dge_genotype, 
          file = "results/differential-gene-expression/dge_genotype.csv",
          row.names = F)

# significant DEG per cluster
dge_sig <- dge_genotype[dge_genotype$p_val_adj <=0.05,]
dge_sig$fold <- abs(dge_sig$avg_logFC)
dge_sig <- dge_sig[dge_sig$fold >= 0.25,]
dge_sig$regulation <- factor(dge_sig$regulation, levels = c("Up", "Down"))
dge_sig$cluster <- factor(dge_sig$cluster, levels = c("SMC-Mod1",
                                                      "SMC-Mod2",
                                                      "SMC-Mod3",
                                                      "SMC-Mod4",
                                                      "SMC-Con",
                                                      "EC",
                                                      "MAC"))

# plot
p <- ggplot(dge_sig, aes(x=(cluster))) +
  geom_bar(aes(fill = regulation),
           position = position_stack(reverse = F)) +
  scale_fill_manual(name = "Regulation in KO",
                    values = c("#D75438", "#4878CD")) +
  xlab("") +
  ylim(0,120) +
  RotatedAxis() +
  ylab("Number of differentially expressed genes") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_blank())

# export
pdf(file = "results/differential-gene-expression/dge_per_cluster.pdf",
    height = 5,
    width = 5,
    useDingbats = F)
print(p)
dev.off()












## DGE KO SMC-Mods v WT SMC-Mods----
dge_mod_genotype <- FindMarkers(obj, 
                                subset.ident = c("SMC-Mod1",
                                                 "SMC-Mod2",
                                                 "SMC-Mod3",
                                                 "SMC-Mod4"),
                                group.by = "genotype",
                                ident.1 = "KO",
                                assay = "RNA",
                                slot = "data",
                                logfc.threshold = 0,
                                verbose = T)
dge_mod_genotype$gene <- rownames(dge_mod_genotype)
dge_mod_genotype$regulation <- ifelse(dge_mod_genotype$avg_logFC > 0, "Up", "Down")
write.csv(dge_mod_genotype,
           file = "results/differential-gene-expression/dge_mod_genotype.csv",
           row.names = F)

## DGE Contractile v. Modulated----
dge_mod_con <- FindMarkers(obj, 
                           ident.1 = c("SMC-Mod1",
                                       "SMC-Mod2",
                                       "SMC-Mod3",
                                       "SMC-Mod4"),
                           ident.2 = "SMC-Con",
                           assay = "RNA",
                           slot = "data",
                           logfc.threshold = 0,
                           verbose = T)
dge_mod_con$gene <- rownames(dge_mod_con)
dge_mod_con$regulation <- ifelse(dge_mod_con$avg_logFC > 0, "Up", "Down")
write.csv(dge_mod_con,
          file = "results/differential-gene-expression/dge_mod_con.csv",
          row.names = F)
