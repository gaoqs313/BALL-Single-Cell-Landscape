library(Seurat)

## Load meta table with sample ID and subtype
meta <- read.table("meta.txt", sep = "\t")

## Merge samples (data normalized in the same way, so set merge.data = TRUE)
x <- readRDS(paste0(meta[1,1], "/obj.rds"))
y <- vector()
for(i in 2:nrow(meta))
{
        assign(meta[i,1], readRDS(paste0(meta[i,1], "/obj.rds")))
        y <- c(y, get(meta[i,1]))
}
merged <- merge(x = x, y = y, add.cell.ids = meta[,1], project = "BALL", merge.data = TRUE)

## Perform cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged <- CellCycleScoring(merged, g2m.features = g2m.genes, s.features = s.genes)

## Normalize
merged <- NormalizeData(merged)

## Find variable
merged <- FindVariableFeatures(merged, selection.method = "vst")

## Scale
vars_to_regress <- c("CC.Difference", "percent.mito", "nCount_RNA", "nFeature_RNA")
merged <- ScaleData(merged, vars.to.regress = vars_to_regress)

## Run PCA
merged <- RunPCA(merged, npcs = 200)

## Run UMAP/tSNE
component <- 120
merged <- RunUMAP(merged, dims = 1:component)

## Find clusters
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:component)
merged <- FindClusters(merged, resolution = 2)

## Specify subtype color
subtypeCol=c()
{
  subtypeCol["DUX4"]='grey40'
  subtypeCol["ETV6-RUNX1-like"]="deeppink"
  subtypeCol["Hyperdiploid"]="#3E9F32"
  subtypeCol["iAMP21"]="lightslateblue"
  subtypeCol["KMT2A"]="#1F78B5"
  subtypeCol["Low-hypodiploid"]="#1E90FF"
  subtypeCol["MEF2D"]="#66C2A6"
  subtypeCol["Near-haploid"]='#B8B828'
  subtypeCol["Other"]='grey75'
  subtypeCol["BCR::ABL1"]="magenta3"
  subtypeCol["BCR::ABL1-like"]="brown"
  subtypeCol["PAX5alt"]="#FFA620"
  subtypeCol["TCF3-PBX1"]="darkgoldenrod4"
  subtypeCol["ZNF384"]="#A8DD00"
}

## Specify cell type color
celltypeCol=c()
{
  celltypeCol["B"]='#4DAF4A'    
  celltypeCol["Blast"]='#377EB8' 
  celltypeCol["DC"]='#FF7F00'    
  celltypeCol["Erythroid"]='#E41A1C'  
  celltypeCol["HSC/MPP and pro."]='#00CDCD' 
  celltypeCol["Monocyte"]='#984EA3' 
  celltypeCol["Plasmablast"]='#EBBC2E' 
  celltypeCol["T_NK"]='#A65628'    
}

## Plot UMAP by cell type, subtype, tissue, or sample
pdf("Final_UMAP_Merged_CellType.pdf", width = 8.7)
DimPlot(merged, reduction = "umap", label = F, pt.size = 0.1, group.by = "CellType", cols = celltypeCol, raster = T) + ggtitle("")
dev.off()

pdf("Final_UMAP_Merged_Subtype.pdf", width = 8.8)
DimPlot(merged, reduction = "umap", label = F, pt.size = 0.1, group.by = "Subtype", cols=subtypeCol, raster = T) + ggtitle("")
dev.off()

pdf("Final_UMAP_Merged_Tissue.pdf", width = 8.7)
DimPlot(merged, reduction = "umap", label = F, pt.size = 0.1, group.by = "Tissue", raster = T) + ggtitle("")
dev.off()

pdf("Final_UMAP_Merged_Sample.pdf", width = 16.5)
DimPlot(merged, reduction = "umap", label = F, pt.size = 0.1, group.by = "Sample", raster = T) + ggtitle("")
dev.off()
