library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(hdf5r)

#Setting working directory
setwd("~/T10")

#Importing data from files
rawData1 <- Read10X_h5(filename ="data/1_JEPDX1_6/sample_filtered_feature_bc_matrix.h5")
rawData2 <- Read10X_h5(filename ="data/2_CHPDX1_5/sample_filtered_feature_bc_matrix.h5")
rawData3 <- Read10X_h5(filename ="data/3_JEPDX2_5/sample_filtered_feature_bc_matrix.h5")
rawData4 <- Read10X_h5(filename ="data/4_CHPDX2_5/sample_filtered_feature_bc_matrix.h5")
rawData5 <- Read10X_h5(filename ="data/5_CHPDX2_7/sample_filtered_feature_bc_matrix.h5")
rawData6 <- Read10X_h5(filename ="data/6_JEPDX3_5/sample_filtered_feature_bc_matrix.h5")
rawData7 <- Read10X_h5(filename ="data/7_CHPDX3_3/sample_filtered_feature_bc_matrix.h5")
rawData8 <- Read10X_h5(filename ="data/8_CHPDX3_5/sample_filtered_feature_bc_matrix.h5")

#Creating SEURAT objects
serData1 <- CreateSeuratObject(rawData1, project = "Mouse1")
serData2 <- CreateSeuratObject(rawData2, project = "Mouse2")
serData3 <- CreateSeuratObject(rawData3, project = "Mouse3")
serData4 <- CreateSeuratObject(rawData4, project = "Mouse4")
serData5 <- CreateSeuratObject(rawData5, project = "Mouse5")
serData6 <- CreateSeuratObject(rawData6, project = "Mouse6")
serData7 <- CreateSeuratObject(rawData7, project = "Mouse7")
serData8 <- CreateSeuratObject(rawData8, project = "Mouse8")

#Adding patient information
serData1$type <- "patientSample1"
serData2$type <- "patientSample1"
serData3$type <- "patientSample2"
serData4$type <- "patientSample2"
serData5$type <- "patientSample2"
serData6$type <- "patientSample3"
serData7$type <- "patientSample3"
serData8$type <- "patientSample3"

#Merging into 1 object
allData <- merge(serData1, c(serData2, serData3, serData4, serData5, serData6, serData7, serData8), add.cell.ids = c("mouse1", "mouse2", "mouse3", "mouse4", "mouse5", "mouse6", "mouse7", "mouse8"))
allData <- JoinLayers(allData)

#Remove unecessary objects
rm(rawData1, rawData2, rawData3, rawData4, rawData5, rawData6, rawData7, rawData8, serData1, serData2, serData3, serData4, serData5, serData6, serData7, serData8)
gc()

#Showing Mitochondrial, Ribosomal, Hemoglobin, Platelate markers
allData <- PercentageFeatureSet(allData, "^MT-", col.name = "percent_mito")
allData <- PercentageFeatureSet(allData, "^RP[SL]", col.name = "percent_ribo")
allData <- PercentageFeatureSet(allData, "^HB[^(P|E|S)]", col.name = "percent_hb")
allData <- PercentageFeatureSet(allData, "PECAM1|PF4", col.name = "percent_plat")

#Filter out cells and genes with low expression
selectCell <- WhichCells(allData, expression = nFeature_RNA > 100)
selectFeat <- rownames(allData)[Matrix::rowSums(allData[["RNA"]]$counts) > 2]
data.filt <- subset(allData, features = selectFeat, cells = selectCell)
dim(data.filt)
table(data.filt$orig.ident)

#Filter out cells with high mitocondrial and low ribosomal expression
data.filt <- subset(data.filt, percent_mito < 30 & percent_ribo > 5)
dim(data.filt)
table(data.filt$orig.ident)

# Filter Mitocondrial, Ribosomal and Hemoglobin genes
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
data.filt <- data.filt[ ! grepl("^RP[SL]", rownames(data.filt)), ]
data.filt <- data.filt[!grepl("^HB[^(P|E|S)]", rownames(data.filt)), ]
dim(data.filt)

data.filt[["RNA"]] <- split(data.filt[["RNA"]], f = data.filt$orig.ident)

data.norm <- NormalizeData(data.filt)
data.foundvar <- FindVariableFeatures(data.norm)
data.scaled <- ScaleData(data.foundvar)
data.pca <- RunPCA(data.scaled, verbose = FALSE)

#Integrate
data.integ <- IntegrateLayers(object = data.pca, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",
                             verbose = FALSE)
data.integ[["RNA"]] <- JoinLayers(data.integ[["RNA"]])


#Plot global variability
data.umap <- RunUMAP(data.integ, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony")
DimPlot(data.umap, reduction = "umap_harmony", group.by = "orig.ident") + ggtitle("Harmony UMAP")

data.umapNonInt <- RunUMAP(data.pca, dims = 1:30, reduction = "pca")
DimPlot(data.umapNonInt, reduction = "umap", group.by = "orig.ident") + ggtitle("Harmony UMAP NonIntegrated")

#Plot with in group variability
data.tsne <- RunTSNE(
  data.integ,
  reduction = "harmony", dims = 1:15,
  perplexity = 30,
  max_iter = 1000,
  theta = 0.5,
  eta = 200,
  num_threads = 0
)
DimPlot(data.tsne, reduction = "tsne", group.by = "orig.ident")

#Plot with in group variability
data.tsneNonInt <- RunTSNE(
  data.pca,
  reduction = "pca", dims = 1:15,
  perplexity = 30,
  max_iter = 1000,
  theta = 0.5,
  eta = 200,
  num_threads = 0
)
DimPlot(data.tsneNonInt, reduction = "tsne", group.by = "orig.ident")

#Summorize based on mouse
wrap_plots(
  DimPlot(data.umapNonInt, reduction = "umap", group.by = "orig.ident") + ggtitle("NonIntegrated UMAP"),
  DimPlot(data.umap, reduction = "umap_harmony", group.by = "orig.ident") + ggtitle("Harmony UMAP"),
  DimPlot(data.tsneNonInt, reduction = "tsne", group.by = "orig.ident") + ggtitle("NonIntegrated tSNE"),
  DimPlot(data.tsne, reduction = "tsne", group.by = "orig.ident") + ggtitle("Harmony tSNE"),
  ncol = 2
) + plot_layout(guides = "collect")

wrap_plots(
  DimPlot(data.umapNonInt, reduction = "umap", group.by = "type") + ggtitle("NonIntegrated UMAP"),
  DimPlot(data.umap, reduction = "umap_harmony", group.by = "type") + ggtitle("Harmony UMAP"),
  DimPlot(data.tsneNonInt, reduction = "tsne", group.by = "type") + ggtitle("NonIntegrated tSNE"),
  DimPlot(data.tsne, reduction = "tsne", group.by = "type") + ggtitle("Harmony tSNE"),
  ncol = 2
) + plot_layout(guides = "collect")

data.umap <- FindNeighbors(data.umap, dims = 1:20)
data.umap <- FindClusters(data.umap, resolution = 0.5)

data.tsne <- FindNeighbors(data.tsne, dims = 1:20)
data.tsne <- FindClusters(data.tsne, resolution = 0.5)

wrap_plots(
  DimPlot(data.umap, reduction = "umap_harmony") + ggtitle("UMAP clustered"),
  DimPlot(data.umap, reduction = "umap_harmony", group.by = "type") + ggtitle("UMAP patientSample"),
  DimPlot(data.tsne, reduction = "tsne") + ggtitle("tSNE clustered"),
  DimPlot(data.tsne, reduction = "tsne", group.by = "type") + ggtitle("tSNE patientSample"),
  ncol = 2
) + plot_layout(guides = "collect")

