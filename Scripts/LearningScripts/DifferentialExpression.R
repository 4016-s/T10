library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(hdf5r)
library(SeuratData)
library(ggrepel)

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

#Filter out cells with high mitocondrial and low ribosomal expression
data.filt <- subset(data.filt, percent_mito < 30 & percent_ribo > 5)

# Filter Mitocondrial, Ribosomal and Hemoglobin genes
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
data.filt <- data.filt[ ! grepl("^RP[SL]", rownames(data.filt)), ]
data.filt <- data.filt[!grepl("^HB[^(P|E|S)]", rownames(data.filt)), ]

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

#Clustering
data.umap <- FindNeighbors(data.umap, dims = 1:30)
data.umap <- FindClusters(data.umap, resolution = 0.5)

DimPlot(data.umap, reduction = "umap_harmony") + ggtitle("UMAP clustered")

#pseudobulk
data.pseudo <- AggregateExpression(data.umap, assays = "RNA", return.seurat = T, group.by = c("type", "orig.ident"))

Idents(data.pseudo) <- "type"

data.de23 <- FindMarkers(object = data.pseudo, 
                            ident.1 = "patientSample2", 
                            ident.2 = "patientSample3",
                            test.use = "DESeq2")
head(data.de23, n = 15)
#significant genes
data.de23$gene <- rownames(data.de23)
signGenes <- data.de23$gene[which(data.de23$p_val_adj < 0.05)]
print(length(signGenes))

data.umap$mouse_type <- paste0(data.umap$orig.ident, "-", data.umap$type)

# generate violin plot 
Idents(data.umap) <- "type"
VlnPlot(data.umap, features = signGenes[1:2], idents = c("patientSample2", "patientSample3"), group.by = "type") 
VlnPlot(data.umap, features = signGenes[1:2], idents = c("patientSample2", "patientSample3"), group.by = "mouse_type", ncol = 2)

FeaturePlot(data.umap, features = signGenes[1:9])
DoHeatmap(data.umap, features = signGenes[1:100], group.by = "mouse_type")

#Volcanoplot
volcplot.de23 <- data.de23 %>%
  mutate(gene = rownames(data.de23),
         significance = case_when(
           avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Upregulated",
           avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Downregulated",
           TRUE ~ "Not significant"
         ))
top_genes <- volcplot.de23 %>%
  arrange(p_val_adj) %>%
  slice(1:20)

ggplot(volcplot.de23, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Downregulated" = "#00AFBB",
                                "Not significant" = "grey",
                                "Upregulated" = "#bb0c00")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_point(size = 1) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 2, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "Pat 2 vs 3: Differentially Expressed Genes",
       x = "Average log2 Fold Change",
       y = "-log10 Adjusted P-value")
