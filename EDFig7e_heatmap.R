
library(ComplexHeatmap)
library(circlize)

## Load top genes
gtop <- read.table("program.Zscore.top30gene.txt", header = T)

## Extract top 30 genes in C and D
topD <- gtop[, 1]
topE <- gtop[, 2]

## Load samples
sampleInfor=read.table("/home/qgao/zdrive/Projects/scRNA_w_Ilaria/Reference/BulkRNA/1.Count/old/input.STD4tSNE2048S.info.txt", header = T, as.is = T, sep = "\t", na.strings = "notUsedString")

## Load DUX4 subtype (only stranded total RNA is used)
dux4 <- read.table("DUX4_sub_bulk.txt", row.names = 1)
dux4$library <- sampleInfor$library[match(rownames(dux4), sampleInfor$sample)]

## Load rlog
rlog <- read.table("/home/qgao/zdrive/Projects/scRNA_w_Ilaria/Reference/BulkRNA/3.rlog/old/rlog_gene_2048S.txt")

## Load program score
score <- read.table("ssGSEA_program_score_by_top30_signature_genes_bulk.txt", header = T, sep = '\t')

## Subset
score_subset <- score[score$Sample %in% rownames(dux4), c(1,2,3)]
score_subset$diff <- score_subset$A - score_subset$B
score_subset$diff <- scale(score_subset$A) - scale(score_subset$B)
score_subset <- score_subset[order(-score_subset$diff), ]

## Subset rlog
df <- rlog[match(c(topD, topE), rownames(rlog)), rownames(dux4)]

## Sort
df <- df[, match(score_subset$Sample, colnames(df))]

## Scale by genes
scaled_rlog = t(scale(t(df)))

## Prepare annotation
anno <- dux4[match(colnames(scaled_rlog), rownames(dux4)),,drop=F]
ha = HeatmapAnnotation(Score = anno_barplot(score_subset$diff, bar_width=1, ylim=c(-4, 4), gp=gpar(fill=ifelse(score_subset$diff<0, "#B2182B", "#2166AC"), col=ifelse(score_subset$diff<0, "#B2182B", "#2166AC"))), Type=anno[,1], Library=anno[,2], col=list(Type=c("D2"="firebrick", "D1"="dodgerblue"), Library=c('Stranded_total_RNA_PE100bp'='#66C2A5', 'Stranded_total_RNA_PE125bp'='#8DA0CB', 'Stranded_total_RNA_PE75bp'='#A6D854', 'Unstranded_mRNA_PE100bp'='#FC8D62', 'Unstranded_mRNA_PE75bp'='#E78AC3')), annotation_height = unit(c(20,4,4), "mm"), annotation_name_side='left', annotation_legend_param = list(direction = "horizontal"), annotation_label=c(paste("ssGSEA", "Score", "Difference", sep="\n"), "Type", "Library"))

## Output
df_out <- cbind(score_subset, Type=anno[,1])
write.table(df_out, "DUX4_programe_score_to_check_with_outcome.txt", row.names=F, col.names=T, sep="\t", quote=F)

## Specify heatmap color
myColor <- colorRamp2(seq(-2,2,0.04), colorRampPalette(colors = c("blue4","blue","white","red","red4") )(101))

## Change rownames
rownames(scaled_rlog) <- rep('', 60)
rownames(scaled_rlog)[c(15,45)] <- c("Stemness", "Inflammation")

## plot
pdf("DUX4_expression_heatmap_color.pdf", width=5.5, height=5, useDingbats=F)
ht <- Heatmap(scaled_rlog, col=myColor, cluster_rows=F, cluster_columns=F, show_row_names=T, row_names_side='left', show_column_names=F, show_row_dend=F, show_column_dend=F, show_heatmap_legend=T, name='Expression', heatmap_legend_param = list(direction = "horizontal"), border=T, top_annotation=ha)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend=T, padding = unit(c(8, 2, 2, 2), "mm"))
decorate_heatmap_body("Expression", {
    grid.lines(c(0, 1), c(1/2, 1/2), gp = gpar( lwd = 1, col = 1))
})
dev.off()


