#ANALYSIS FOR AIM+ CD8
#Load libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/scRNA-seq_AIM_649/CD8_GEX.MPS12349973-D12.sorted.2155.1.10x_outputs")

#Create Seurat objects
counts <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
seurat_CD8 <- CreateSeuratObject(counts, project = "CD8")
rm(counts)

#Quality control
seurat_CD8[["percent.mt"]] <- PercentageFeatureSet(seurat_CD8, "^MT-")
VlnPlot(seurat_CD8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_CD8 <- subset(seurat_CD8, 
                     subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
seurat_CD8

#Normalization
seurat_CD8 <- NormalizeData(seurat_CD8)

#Feature selection
seurat_CD8 <- FindVariableFeatures(seurat_CD8)

top10 <- head(VariableFeatures(seurat_CD8), 10)
plot1 <- VariableFeaturePlot(seurat_CD8)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

#Scaling
seurat_CD8 <- ScaleData(seurat_CD8)

#PCA
seurat_CD8 <- RunPCA(seurat_CD8, npcs = 50)

#UMAP
seurat_CD8 <- RunUMAP(seurat_CD8, reduction = "pca", dims = 1:20)
DimPlot(seurat_CD8, reduction = "umap")

#Clustering
seurat_CD8 <- FindNeighbors(seurat_CD8, dims = 1:20) 
seurat_CD8 <- FindClusters(seurat_CD8, resolution = 0.8)
table(seurat_CD8@meta.data$seurat_clusters)

for(res in c(0.2, 0.5, 1, 1.5, 2)) {
  seurat_CD8<-FindClusters(seurat_CD8, resolution=res)
  print(table(seurat_CD8@meta.data$seurat_clusters))
}

#Finding optimal clusters with Silhouette Index
require(cluster) # Contains silhouette function
dist.matrix <- dist(x = Embeddings(object = seurat_CD8[["pca"]])[, 1:20])
sil <- silhouette(as.numeric(as.character(seurat_CD8@meta.data$RNA_snn_res.0.8)), dist=dist.matrix)
head(sil)
mean(sil[,3])

cluster_sil_scores <- aggregate(sil[,3], by=list(sil[,1]), mean)
cluster_sil_scores

clusterings <- grepl("RNA_snn_res",colnames(seurat_CD8@meta.data))
clusterings <- colnames(seurat_CD8@meta.data)[clusterings]
for (c in clusterings) {
  require(cluster) 
  sil <- silhouette(as.numeric(as.character(seurat_CD8@meta.data[,c])), dist=dist.matrix)
  print(mean(sil[,3]))
}
clusterings # To check the order of the resolutions
seurat_CD8 <- FindClusters(seurat_CD8, resolution = 0.8) # Make this resolution default

saveRDS(seurat_CD8, "./seurat_CD8.rds")

#Visualize clustering
DimPlot(seurat_CD8, reduction = "umap", label = T, pt.size = 1.5, label.size = 5)

#Cell cluster annotation
##Feature plot of canonical markers
selected_markers_CD8 <- c("SELL", "TCF7", # Naive
                 "IL2RA", "IFNG", # Effector
                 "IL7R", "CCR7", # Memory
                 "TOX", "PDCD1") # Exhausted
FeaturePlot(seurat_CD8, features = selected_markers_CD8, 
            ncol = 2, reduction = "umap", cols = c("seashell2", "brown2")) & NoAxes() & NoLegend()

##Dot plot of canonical markers
markers_CD8 <- c("SELL", "CCR7", "LEF1", "TCF7",
                 "KLRG1", "IL2RA", "TBX21", "IFNG", "TNF", "GZMB", "PRF1", "CCL3", "CCL4",
                 "IL7R", "EOMES", "GZMK",
                 "CX3CR1", "FGFBP2",
                 "TOX", "PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "ENTPD1",
                 "SLC4A10","TRAV1-2", "ZBTB16", "KLRB1")
library(rlist)
markers_CD8_reverse <- list.reverse(markers_CD8)
DotPlot(seurat_CD8, features = markers_CD8_reverse, cluster.idents = T, 
        cols = c("white", "red")) + RotatedAxis() + coord_flip() + labs(x = "", y = "")

##Heatmap of marker genes per cluster
all_markers_CD8 <- FindAllMarkers(seurat_CD8, only.pos = TRUE)
all_markers_CD8 <- all_markers_CD8 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) 

saveRDS(all_markers_CD8, "./all_markers_CD8.rds")

table(all_markers_CD8$cluster)
View(all_markers_CD8)
library(writexl)
write_xlsx(all_markers_CD8, "./all_markers_CD8.xlsx")

top10 <- all_markers_CD8 %>%  
  slice_head(n = 10) %>%
  ungroup()

DoHeatmap(seurat_CD8,features = top10$gene, size = 3)

##Rename clusters
new_identities <- c("EM I", "Eff I", "Eff III", "TEMRA", "Exh", 
                    "Eff II", "Exh-Pre", "Exh-Term",
                    "ISG", "Eff IV", "EM II")
names(new_identities) <- levels(seurat_CD8)
seurat_CD8 <- RenameIdents(seurat_CD8, new_identities)

DimPlot(seurat_CD8, reduction = "umap", label = T, pt.size = 1.5, label.size = 5) + NoLegend()
DotPlot(seurat_CD8, features = markers_CD8_reverse, cluster.idents = T, 
        cols = c("white", "red")) + RotatedAxis() + coord_flip() + labs(x = "", y = "")
DoHeatmap(seurat_CD8,features = top10$gene, size = 0)

#VDJ analysis
library(scRepertoire)
library(ggplot2)
library(Seurat)

setwd("~/Desktop/scRNA-seq_AIM_649/CD8_VDJ.MPS12349777-F10.sorted.2155.1.10x_outputs")

##Loading and processing contig data
s1 <- read.csv("./filtered_contig_annotations.csv")

contig_list <- list(s1)
head(contig_list[[1]])

##Combining contigs into clones
combined.TCR <- combineTCR(contig_list,
                           samples = "CD8",
                           removeNA = FALSE,
                           removeMulti = FALSE,
                           filterMulti = FALSE)

head(combined.TCR[[1]])

##Export the paired clonotypes
exportClones(combined.TCR,
             write.file = TRUE,
             dir = "./",
             file.name = "clones.csv")

##Number of unique clonotypes
clonalQuant(combined.TCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = FALSE)

##Relative distribution of clonotypes by abundance
clonalAbundance(combined.TCR,
                cloneCall = "strict",
                scale = FALSE)

##Length distribution of the CDR3 sequences
clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "both")

##Clonal Space Homeostasis
colorblind_vector <- hcl.colors(n = 8, palette = "inferno", fixup = TRUE)

clonalHomeostasis(combined.TCR, 
                  cloneCall = "strict") + 
  scale_fill_manual(values = rev(colorblind_vector[c(8,6,7,5,1)]))

##Relative usage of genes
plot1 <- vizGenes(combined.TCR,
                  x.axis = "TRAV",
                  y.axis = NULL,
                  plot = "barplot",
                  scale = TRUE) + 
  theme(axis.text=element_text(size=12, face = "bold"))

plot2 <- vizGenes(combined.TCR,
                  x.axis = "TRBV",
                  y.axis = NULL,
                  plot = "barplot",
                  scale = TRUE) + 
  theme(axis.text=element_text(size=12, face = "bold"))

plot1 / plot2

##Proportion of amino acids along the cdr3 sequence
percentAA(combined.TCR,
          chain = "TRA",
          aa.length = 20)

##Combining clones and Seurat object
seur_obj <- seurat_CD8
DimPlot(seur_obj, reduction = "umap", label = T, pt.size = 1.5, label.size = 5) + NoLegend()
table(Idents(seur_obj))
head(seur_obj)

cell.barcodes <- rownames(seur_obj[[]])
cell.barcodes <- paste0(seur_obj$orig.ident, "_", cell.barcodes)
seur_obj <- RenameCells(seur_obj, new.names = cell.barcodes)
head(seur_obj)

seur_obj <- combineExpression(combined.TCR, 
                              seur_obj, 
                              cloneCall="strict",
                              proportion = TRUE)

Seurat::DimPlot(seur_obj, group.by = "cloneSize", pt.size = 1.5) +
  scale_color_manual(values = rev(colorblind_vector[c(5,7,6,8)])) + NoAxes() + labs(title = "")

###############################################################################################
#ANALYSIS FOR AIM+ CD4
#Load libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/scRNA-seq_AIM_649/CD4_GEX.MPS12349973-E12.sorted.2155.1.10x_outputs")

#Create Seurat objects
counts <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
seurat_CD4 <- CreateSeuratObject(counts, project = "CD4")
rm(counts)

#Quality control
seurat_CD4[["percent.mt"]] <- PercentageFeatureSet(seurat_CD4, "^MT-")
VlnPlot(seurat_CD4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_CD4 <- subset(seurat_CD4, 
                     subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
seurat_CD4

#Normalization
seurat_CD4 <- NormalizeData(seurat_CD4)

#Feature selection
seurat_CD4 <- FindVariableFeatures(seurat_CD4)

top10 <- head(VariableFeatures(seurat_CD4), 10)
plot1 <- VariableFeaturePlot(seurat_CD4)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

#Scaling
seurat_CD4 <- ScaleData(seurat_CD4)

#PCA
seurat_CD4 <- RunPCA(seurat_CD4, npcs = 50)

#UMAP
seurat_CD4 <- RunUMAP(seurat_CD4, reduction = "pca", dims = 1:20)
DimPlot(seurat_CD4, reduction = "umap")

#Clustering
seurat_CD4 <- FindNeighbors(seurat_CD4, dims = 1:20) 
seurat_CD4 <- FindClusters(seurat_CD4, resolution = 0.8)
table(seurat_CD4@meta.data$seurat_clusters)

for(res in c(0.2, 0.5, 1, 1.5, 2)) {
  seurat_CD4<-FindClusters(seurat_CD4, resolution=res)
  print(table(seurat_CD4@meta.data$seurat_clusters))
}

#Finding optimal clusters with Silhouette Index
require(cluster) # Contains silhouette function
dist.matrix <- dist(x = Embeddings(object = seurat_CD4[["pca"]])[, 1:20])
sil <- silhouette(as.numeric(as.character(seurat_CD4@meta.data$RNA_snn_res.0.8)), dist=dist.matrix)
head(sil)
mean(sil[,3])

cluster_sil_scores <- aggregate(sil[,3], by=list(sil[,1]), mean)
cluster_sil_scores

clusterings <- grepl("RNA_snn_res",colnames(seurat_CD4@meta.data))
clusterings <- colnames(seurat_CD4@meta.data)[clusterings]
for (c in clusterings) {
  require(cluster) 
  sil <- silhouette(as.numeric(as.character(seurat_CD4@meta.data[,c])), dist=dist.matrix)
  print(mean(sil[,3]))
}
clusterings # To check the order of the resolutions
seurat_CD4 <- FindClusters(seurat_CD4, resolution = 0.8) # Make this resolution default

saveRDS(seurat_CD4, "./seurat_CD4.rds")

#Visualize clustering
DimPlot(seurat_CD4, reduction = "umap", label = T, pt.size = 1.5, label.size = 5)

#Cell cluster annotation
##Feature plot of canonical markers
selected_markers_CD4 <- c("TBX21", "CXCR3", # Th1
                          "GATA3", "IL13", # Th2
                          "RORC", "IL17A", # Th17
                          "PDCD1", "CXCR5", # Tfh
                          "FOXP3", "CTLA4") # Treg
FeaturePlot(seurat_CD4, features = selected_markers_CD4, 
            ncol = 2, reduction = "umap", cols = c("seashell2", "mediumpurple4")) & NoAxes() & NoLegend()

##Dot plot of canonical markers
markers_CD4 <- c("TBX21", "CXCR3", "CCR5", "IFNG",
                 "GATA3", "STAT5A", "CCR3", "IL4", "IL5", "IL13",
                 "RORC", "CCR4", "CCR6", "IL17A",
                 "BCL6", "PDCD1", "CXCR5", "IL21",
                 "FOXP3", "IL10", "CTLA4", "IL2RA",
                 "SELL", "IL7R", "CCR7", "TCF7")

library(rlist)
markers_CD4_reverse <- list.reverse(markers_CD4)
DotPlot(seurat_CD4, features = markers_CD4_reverse, cluster.idents = T, 
        cols = c("white", "mediumpurple4")) + RotatedAxis() + coord_flip() + labs(x = "", y = "")

##Heatmap of marker genes per cluster
all_markers_CD4 <- FindAllMarkers(seurat_CD4, only.pos = TRUE)
all_markers_CD4 <- all_markers_CD4 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) 

saveRDS(all_markers_CD4, "./all_markers_CD4.rds")

table(all_markers_CD4$cluster)
View(all_markers_CD4)
library(writexl)
write_xlsx(all_markers_CD4, "./all_markers_CD4.xlsx")

top10 <- all_markers_CD4 %>%  
  slice_head(n = 10) %>%
  ungroup()

DoHeatmap(seurat_CD4,features = top10$gene, size = 3)

##Rename clusters
new_identities <- c("Tfh I", "Th1", "Treg/Th17", "Tfh II", "Th2", "Naive", "Treg")
names(new_identities) <- levels(seurat_CD4)
seurat_CD4 <- RenameIdents(seurat_CD4, new_identities)

DimPlot(seurat_CD4, reduction = "umap", label = T, pt.size = 1.5, label.size = 5) + NoLegend()
DotPlot(seurat_CD4, features = markers_CD4_reverse, cluster.idents = T, 
        cols = c("white", "mediumpurple4")) + RotatedAxis() + coord_flip() + labs(x = "", y = "")
DoHeatmap(seurat_CD4,features = top10$gene, size = 0)

#VDJ analysis
library(scRepertoire)
library(ggplot2)
library(Seurat)

setwd("~/Desktop/scRNA-seq_AIM_649/CD4_VDJ.MPS12349777-H10.sorted.2155.1.10x_outputs")

##Loading and processing contig data
s1 <- read.csv("./filtered_contig_annotations.csv")

contig_list <- list(s1)
head(contig_list[[1]])

##Combining contigs into clones
combined.TCR <- combineTCR(contig_list,
                           samples = "CD4",
                           removeNA = FALSE,
                           removeMulti = FALSE,
                           filterMulti = FALSE)

head(combined.TCR[[1]])

##Export the paired clonotypes
exportClones(combined.TCR,
             write.file = TRUE,
             dir = "./",
             file.name = "clones.csv")

##Number of unique clonotypes
clonalQuant(combined.TCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = FALSE)

##Relative distribution of clonotypes by abundance
clonalAbundance(combined.TCR,
                cloneCall = "strict",
                scale = FALSE)

##Length distribution of the CDR3 sequences
clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "both")

##Clonal Space Homeostasis
colorblind_vector <- hcl.colors(n = 8, palette = "inferno", fixup = TRUE)

clonalHomeostasis(combined.TCR, 
                  cloneCall = "strict") + 
  scale_fill_manual(values = rev(colorblind_vector[c(8,6,7,5,1)]))

##Relative usage of genes
plot1 <- vizGenes(combined.TCR,
                  x.axis = "TRAV",
                  y.axis = NULL,
                  plot = "barplot",
                  scale = TRUE) + 
  theme(axis.text=element_text(size=12, face = "bold"))

plot2 <- vizGenes(combined.TCR,
                  x.axis = "TRBV",
                  y.axis = NULL,
                  plot = "barplot",
                  scale = TRUE) + 
  theme(axis.text=element_text(size=12, face = "bold"))

plot1 / plot2

##Proportion of amino acids along the cdr3 sequence
percentAA(combined.TCR,
          chain = "TRA",
          aa.length = 20)

##Combining clones and Seurat object
seur_obj <- seurat_CD4
DimPlot(seur_obj, reduction = "umap", label = T, pt.size = 1.5, label.size = 5) + NoLegend()
table(Idents(seur_obj))
head(seur_obj)

cell.barcodes <- rownames(seur_obj[[]])
cell.barcodes <- paste0(seur_obj$orig.ident, "_", cell.barcodes)
seur_obj <- RenameCells(seur_obj, new.names = cell.barcodes)
head(seur_obj)

seur_obj <- combineExpression(combined.TCR, 
                              seur_obj, 
                              cloneCall="strict",
                              proportion = TRUE)

Seurat::DimPlot(seur_obj, group.by = "cloneSize", pt.size = 1.5) +
  scale_color_manual(values = rev(colorblind_vector[c(5,7)])) + NoAxes() + labs(title = "")
