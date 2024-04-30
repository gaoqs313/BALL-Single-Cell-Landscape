library(circlize)
library(ComplexHeatmap)

## Load meta program score
meta_score <- read.table("meta_score.txt", check.names=F)

## Calculate correlation
nmf_cor <- cor(meta_score)
nmf_cor_dist <- as.dist(1 - nmf_cor)
hclust <- hclust(nmf_cor_dist, method = "ward.D")
dend = as.dendrogram(hclust)
dend <- reorder(dend, wts = rep(1:5,11), agglo.FUN = "min")

## Specify subtype colors
subtypeCol=c()
{
  subtypeCol["DUX4"]='grey40'
  subtypeCol["ETV6::RUNX1-like"]="deeppink"
  subtypeCol["Hyperdiploid"]="#3E9F32"
  subtypeCol["KMT2A"]="#1F78B5"
  subtypeCol["MEF2D"]="#66C2A6"
  subtypeCol["Near-haploid"]='#B8B828'
  subtypeCol["PAX5alt"]='#FFA620'
  subtypeCol["BCR::ABL1"]="magenta3"
  subtypeCol["BCR::ABL1-like"]="brown"
  subtypeCol["TCF3::PBX1"]="darkgoldenrod4"
  subtypeCol["ZNF384"]="#A8DD00"
}

## Specify meta program colors
programCol=c()
{
  programCol['A'] = '#1B9E77'
  programCol['B'] = '#D95F02'
  programCol['C'] = '#7570B3'
  programCol['D'] = '#E7298A'
  programCol['E'] = '#E6AB02'
}

## Prepare annotation
top_annotation = HeatmapAnnotation(Subtype = rep(names(subtypeCol), each = 5),
                            Program = rep(c('A','B','C','D','E'), 11),
			    col=list(Subtype = subtypeCol, Program = programCol),
                            annotation_name_gp = gpar(fontsize = 7))
row_annotation = rowAnnotation(Subtype = rep(names(subtypeCol), each = 5),
                            Program = rep(c('A','B','C','D','E'), 11),
                            col=list(Subtype = subtypeCol, Program = programCol),
                            annotation_name_gp = gpar(fontsize = 7))

## Prepare heatmap color
myColor <- colorRamp2(seq(-1,1,0.02), colorRampPalette(colors = c("#083160", "#2668AA", "#4794C1", "#94C5DD", "#D2E5EF", "#F7F7F7", "#FCDBC8", "#F2A585", "#D46151", "#B01B2F", "#660220"))(101))

## Plot heatmap
pdf("BALL_NMF_correlation_heatmap.pdf", width = 8.4)
ht <- Heatmap(as.matrix(nmf_cor), cluster_columns=dend, col=myColor, cluster_rows=dend, show_row_dend=T, show_column_dend=T, show_row_names=F, show_column_names=F, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6), name="Correlation", top_annotation=top_annotation, left_annotation=row_annotation)
draw(ht, merge_legend=T)
dev.off()
