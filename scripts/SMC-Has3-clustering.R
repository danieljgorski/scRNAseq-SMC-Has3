## Load libraries----
library(Seurat) # v3.1.5

## Read in 10x outputs and create Seurat objects----
lesion_wt <- Read10X(data.dir = "data/WT_filtered_feature_bc_matrix")
lesion_ko <- Read10X(data.dir = "data/KO_filtered_feature_bc_matrix")
lesion_wt <- CreateSeuratObject(counts = lesion_wt,
                                project = "SMC_Has3_WT_lesion",
                                min.cells = 3,
                                min.features = 200)
lesion_ko <- CreateSeuratObject(counts = lesion_ko,
                                project = "SMC_Has3_KO_lesion",
                                min.cells = 3,
                                min.features = 200)

## Adding metadata----
lesion_wt <- AddMetaData(lesion_wt,
                         metadata = "WT",
                         col.name = "genotype")
lesion_wt <- AddMetaData(lesion_wt,
                         col.name = "tissue",
                         metadata = "Lesion")
lesion_ko <- AddMetaData(lesion_ko,
                         metadata = "KO",
                         col.name = "genotype")
lesion_ko <- AddMetaData(lesion_ko,
                         col.name = "tissue",
                         metadata = "Lesion")

## Merge WT and KO Seurat objects----
obj <- merge(x = lesion_wt,
             y = lesion_ko,
             add.cell.ids = c("lesion_wt", "lesion_ko"))
remove(lesion_wt)
remove(lesion_ko)

## QC filtering----
# First pass clustering showed nFeature subset needs to be >800 or clusters
# are created without unique gene expression, dominated by low nFeature
obj[["percent.mt"]] <- PercentageFeatureSet(object = obj, pattern = "^mt-")
VlnPlot(obj, features = c("percent.mt", "nCount_RNA", "nFeature_RNA"), ncol = 3)
obj[["percent.hemo"]] <- PercentageFeatureSet(object = obj, pattern = "^Hb")
VlnPlot(obj, features = c("percent.hemo"))
sum(obj[["percent.hemo"]])
obj <- subset(x = obj,
              subset = nFeature_RNA > 800 &
                nFeature_RNA < 5000 &
                percent.mt < 10 &
                percent.hemo < 5)

## SCTransform, integration and clustering----
obj_list <- SplitObject(obj, split.by = "orig.ident")
obj_list_sct <- list()
for (i in 1:(length(obj_list))) {
   obj_list_sct[[i]] <- SCTransform(obj_list[[i]],
                                  vars.to.regress = "percent.mt",
                                  return.only.var.genes = F,
                                  verbose = T)
}
obj_features <- SelectIntegrationFeatures(object.list = obj_list_sct,
                                          nfeatures = 3000,
                                          verbose = T)
obj_list_sct <- PrepSCTIntegration(object.list = obj_list_sct,
                               anchor.features = obj_features,
                                  verbose = T)
obj_anchors <- FindIntegrationAnchors(object.list = obj_list_sct,
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

## Add published embeddings for consistency----
load("data/published_embeddings.Rdata")
obj[["published_umap"]] <- CreateDimReducObject(
   embeddings = published_embeddings,
   key = "pubumap_",
   assay = DefaultAssay(obj))

## Save integrated object----
save(obj, file = "results/objects/obj.Rdata")
