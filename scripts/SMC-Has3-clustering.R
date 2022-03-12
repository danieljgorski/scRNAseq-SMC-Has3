## Load libraries----
library(Seurat)

## Read in 10x outputs and create Seurat objects----
lesion_WT <- Read10X(data.dir = "data/WT_sample_filtered_feature_bc_matrix")
lesion_KO <- Read10X(data.dir = "data/KO_sample_filtered_feature_bc_matrix")
lesion_WT <- CreateSeuratObject(counts = lesion_WT, 
                                project= "SMC_Has3_WT_lesion", 
                                min.cells = 3, 
                                min.features = 200)
lesion_KO <- CreateSeuratObject(counts = lesion_KO, 
                                project= "SMC_Has3_KO_lesion", 
                                min.cells = 3, 
                                min.features = 200)

## Adding metadata----
head(lesion_WT@meta.data)
lesion_WT <- AddMetaData(lesion_WT, 
                         metadata = "WT", 
                         col.name = "genotype")
lesion_WT <- AddMetaData(lesion_WT, 
                         metadata = "Lesion", 
                         col.name = "tissue")
lesion_KO <- AddMetaData(lesion_KO, 
                         metadata = "KO", 
                         col.name = "genotype")
lesion_KO <- AddMetaData(lesion_KO, 
                         metadata = "Lesion", 
                         col.name = "tissue")

## Merge WT and KO Seurat objects----
obj <- merge(x = lesion_WT, 
             y = lesion_KO, 
             add.cell.ids = c("Lesion_WT", "Lesion_KO"))
remove(lesion_WT)
remove(lesion_KO)

## QC filtering----
# First pass clustering showed nFeature subset needs to be >800 or clusters 
# are created without unique gene expression, dominated by low nFeature
obj[["percent.mt"]] <- PercentageFeatureSet(object = obj, pattern = "^mt-")
VlnPlot(obj, features =c("percent.mt","nCount_RNA","nFeature_RNA"), ncol=3)
obj[["percent.hemo"]] <- PercentageFeatureSet(object = obj, pattern = "^Hb")
VlnPlot(obj, features =c("percent.hemo"))
sum(obj[["percent.hemo"]])
obj <- subset(x = obj, 
              subset = nFeature_RNA > 800 & 
                nFeature_RNA < 5000 & 
                percent.mt < 10 & 
                percent.hemo < 5)

## SCTransform, integration and clustering----
obj_list <- SplitObject(obj, split.by = "orig.ident")
for (i in 1:length(obj_list)) {
  obj_list[[i]] <- SCTransform(obj_list[[i]], 
                                  vars.to.regress = "percent.mt",
                                  return.only.var.genes = F,
                                  verbose = T)
}
obj_features <- SelectIntegrationFeatures(object.list = obj_list, 
                                          nfeatures = 3000, 
                                          verbose = T)
obj_list <- PrepSCTIntegration(object.list = obj_list, 
                               anchor.features = obj_features, 
                                  verbose = T)
obj_anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                      normalization.method = "SCT", 
                                      anchor.features = obj_features, 
                                      verbose = T)
obj <- IntegrateData(anchorset = obj_anchors, 
                                   normalization.method = "SCT", 
                                   verbose = T)
obj <- RunPCA(object = obj, verbose = T)
ElbowPlot(object = obj)
obj <- FindNeighbors(object = obj, 
                     dims = 1:17, 
                     verbose = T)
obj <- FindClusters(object = obj, 
                    resolution = 0.5, 
                    verbose = T)
obj <- RunUMAP(object = obj, 
               dims = 1:17, 
               verbose = T)
DimPlot(obj, label = T)
DimPlot(obj, split.by = "genotype")
FeaturePlot(obj, features = "eyfp")
obj

## Normalization and scaling of the RNA assay----
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     verbose = T)
all_genes <- rownames(obj)
obj <- ScaleData(obj, 
                 features = all_genes, 
                 vars.to.regress = "percent.mt", 
                 verbose = T)
levels(obj) 

## Saved integrated object----
save(obj, file = "results/objects/obj.Rdata")


