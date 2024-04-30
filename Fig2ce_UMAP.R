library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)

sample <- 'SJE2A063_D1'

## Load sample RDS
load(paste0(sample, ".RData"))

## Plot UMAP color coded by cell type
pdf(paste0(sample, "_celltype.pdf"), width=7.2, useDingbats = F)
DimPlot(obj, reduction = "umap", label = F, repel = T, raster = T, pt.size = 0.5, group.by = "CellType", cols = celltypeCol) + ggtitle("") + theme(axis.text = element_blank(), axis.ticks=element_blank()) + NoLegend()
dev.off()

## Load CNV
infercnv <- read.table(paste0(inferCNV_folder, sample, "/subclusters/infercnv.observation_groupings.txt"))
obj@meta.data$CNV <- infercnv[match(rownames(obj@meta.data), rownames(infercnv)), 1]
obj@meta.data$CNV[obj@meta.data$CNV %in% c(2,5,6,7,8)] <- "C1"
obj@meta.data$CNV[obj@meta.data$CNV %in% c(3,4)] <- "C2"
obj@meta.data$CNV[obj@meta.data$CNV %in% c(1)] <- "C3"
obj@meta.data$CNV[is.na(obj@meta.data$CNV)] <- "Normal"

## Specify color
cnvColor = c()
{
  cnvColor["C1"]="#1B9E77"
  cnvColor["C2"]="#D95F02"
  cnvColor["C3"]="#7570B3"
  cnvColor["Normal"]="grey75"
}

## Plot UMAP color coded by cnv clones
pdf(paste0(sample, "_cnv_cluster.pdf"), width=7.2, useDingbats = F)
DimPlot(obj, reduction = "umap", label = F, repel = F, raster = T, label.size = 5, pt.size = 0.5, group.by = "CNV", cols = cnvColor) + ggtitle("") + theme(axis.text = element_blank(), axis.ticks=element_blank()) + NoLegend()
dev.off()

## Load Fusion
barcodes <- read.table(paste0(sample, "_cells_with_fusion.barcodes.txt"))
obj@meta.data$Fusion <- NA
obj@meta.data$Fusion[rownames(obj@meta.data) %in% barcodes[,1] & obj@meta.data$CellType == "Blast"] <- "Detected"
highlist <- list(Fusion = rownames(obj@meta.data[!is.na(obj@meta.data$Fusion), ]))

## Plot UMAP color coded by fusion status
pdf(paste0(sample, "_fusion.pdf"), width=7.2, useDingbats = F)
DimPlot(obj, reduction = "umap", label = T, repel = T, label.size = 5, pt.size = 0.5, cells.highlight = highlist, cols.highlight = "red",  cols = "grey", group.by = "Fusion") + ggtitle("") + NoLegend() + theme(axis.text = element_blank(), axis.ticks=element_blank()) + NoLegend()
dev.off()
