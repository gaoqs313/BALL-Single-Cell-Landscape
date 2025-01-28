library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

## Load data
clr <- read.csv("Merged_CLR.csv", row.names=1, header=T, check.names=F)
clr <- as.data.frame(t(clr))

## Load meta data
meta <- read.table("meta_data.txt", sep="\t", row.names=1, header=T)
meta <- meta[match(colnames(clr), rownames(meta)), ]

## Prepare annotation
ha <- HeatmapAnnotation(Subtype=meta$Subtype, Time=meta$Time, col=list(Subtype=c("ETV6::RUNX1"="gold2", "ETV6::RUNX1-like"="deeppink"), Time=c("Diagnostic"="grey90", "Post-treatment"="grey30")))

## Prepare color
myColor <- colorRamp2(seq(0,6,0.06), colorRampPalette(colors = brewer.pal(9, "Purples"))(101))

## Plot heatmap
pdf("BD_ER.pdf", height=5, width=6)
ht <- Heatmap(as.matrix(clr), col=myColor, cluster_rows=F, cluster_columns=T, show_row_names=T, show_column_names=T, show_heatmap_legend=T, column_names_gp = gpar(fontsize=8), name='Abundance', top_annotation=ha)
draw(ht, merge_legend=T)
dev.off()
