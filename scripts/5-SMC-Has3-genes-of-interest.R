## Load libraries----
library(Seurat)
library(ggplot2)
library(readr)
## Load Seurat object----
load("results/objects/obj.Rdata")

## Genes of Interest----
GOIs <- c("Acan",
          "Acta2",
          "Adgre1",
          "C3",
          "Ccl2",
          "Cd68",
          "Cdh5",
          "Chad",
          "Chil1",
          "Cnn1",
          "Cytl1",
          "Col1a1",
          "Col2a1",
          "Col3a1",
          "Col5a1",
          "Col5a2",
          "Col9a1",
          "Col15a1",
          "Comp",
          "Dcn",
          "eyfp",
          "Egr1",
          "Fn1",
          "Fmod",
          "Fcgr1",
          "Gdf10",
          "Hapln1",
          "Ibsp", 
          "Itgb1", 
          "Itga8",
          "Jun",
          "Lbp",
          "Lrg1",
          "Lcn2",
          "Lgals3",
          "Lum",
          "Lbp",
          "Ly6a",
          "Mmp3",
          "Mmp13",
          "Myh11",
          "Pecam1",
          "Ptprc",
          "Phactr1",
          "Pdgfrb", 
          "Prg4",
          "S100b",
          "Spp1",
          "Sox9",
          "Runx2",
          "Comp",
          "Chad",
          "Thbs1",
          "Pi16",
          "Tagln",
          "Timp1",
          "Tnc",
          "Tnfrsf11b",
          "Txnip",
          "Trpv4",
          "Vcan",
          "Vcam1",
          "Zfp36l1",
          "Wif1",
          "Itgam",
          "Mrc1"
)

## FeaturePlots----
for (i in GOIs) {
  f <- FeaturePlot(obj, features = i, pt.size = 1, cols = c("grey", "red")) +
    theme(panel.border = element_rect(colour = 'black', size = 1),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.title = element_blank(),
          legend.position = c(0.08, 0.08),
          legend.text = element_text(size = 12),
          legend.direction = "horizontal") + 
    annotate("text", x = -13.5, y = -5, size = 13, label = i, fontface = "italic")
  
  pdf(file = paste0("results/genes-of-interest/FeaturePlot_", i, ".pdf"),
      width = 6,
      height = 6,
      useDingbats = F)
  print(f)
  dev.off()
}

## VlnPlots----
Dimcols <- c("#EDAE3C",
             "#D95030",
             "#FB5938",
             "#D6809D",
             "#508AB2",
             "#6FA7A0",
             "#7572AD")

for (i in GOIs) {
  v <- VlnPlot(obj, features = i, cols = Dimcols, pt.size = 0) +
    RotatedAxis() +
    ggtitle(i) +
    theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 45),
          legend.position = "none",
          axis.title.y = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, size = 22)) +
    ylab("Expression level") 
  
  pdf(file = paste0("results/genes-of-interest/VlnPlot_", i, ".pdf"),
      width = 6,
      height = 6,
      useDingbats = F)
  print(v)
  dev.off()
}


## Gene set signatures----

# read in gene signatures
gene_signatures <- read_csv("data/gene_signatures.csv")
SLRPs <- na.omit(gene_signatures$SLRPs)
FibCol <- na.omit(gene_signatures$`Fibrillar collagen`)

# score the cells
obj <- AddModuleScore(obj, 
                      features = list(SLRPs), 
                      assay = "RNA",
                      name = "SLRP")
obj <- AddModuleScore(obj, 
                      features = list(FibCol), 
                      assay = "RNA",
                      name = "FibCol")

# VlnPlots
pdf(file = "results/genes-of-interest/VlnPlot_SLRP_signature.pdf",
    width = 6,
    height = 6,
    useDingbats = F)
VlnPlot(obj, features = "SLRP1", cols = Dimcols, pt.size = 0) +
  RotatedAxis() +
  ggtitle("SLRP") +
  ylab("Average expression level") +
  theme(plot.title = element_text(hjust = 0.5, size = 35),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, size = 22)) 
dev.off()

pdf(file = "results/genes-of-interest/VlnPlot_FibCol_signature.pdf",
    width = 6,
    height = 6,
    useDingbats = F)
VlnPlot(obj, features = "FibCol1", cols = Dimcols, pt.size = 0) +
  RotatedAxis() +
  ggtitle("Fibrillar collagen") +
  ylab("Average expression level") +
  theme(plot.title = element_text(hjust = 0.5, size = 35),
        legend.position = "none",
        axis.title.y = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, size = 22)) 
dev.off()
