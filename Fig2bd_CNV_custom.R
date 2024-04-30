library(infercnv)
library(ComplexHeatmap)
library(circlize)
library(stringr)

sample <- 'SJE2A063_D1'

## Load normalized data from infercnv folder
df <- read.table(paste0(inferCNV_folder, sample, "/subclusters/infercnv.observations.txt"), check.names=F)

## Reverse to match infercnv order
df <- df[, ncol(df):1]

## Load groupings from infercnv folder
grp <- read.table(paste0(inferCNV_folder, sample, "/subclusters/infercnv.observation_groupings.txt"))
grp <- grp[match(colnames(df), rownames(grp)), 1, drop = F]
grp[, 1] <- factor(grp[, 1], levels = unique(grp[, 1]), ordered = T)

## Load R object from infercnv folder
objs <- readRDS(paste0(inferCNV_folder, sample, "/subclusters/run.final.infercnv_obj"))

## Extract number of genes for each chromosome
ngene <- as.numeric(table(objs@gene_order$chr))
end <- cumsum(as.numeric(table(objs@gene_order$chr)))
start <- c(0, end[1:22])
middle <- start + (end - start) / 2

## Load heatmap thresholds form infercnv folder
brk <- read.table(paste0(inferCNV_folder, sample, "/subclusters/infercnv.heatmap_thresholds.txt"))
brk <- brk[,1]

## Specify heatmap color
myColor <- colorRamp2(brk, colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(brk)))

## Assign infercnv groupings to cnv clones
obj_meta <- grp
colnames(obj_meta) <- "CNV"
obj_meta$CNV <- as.character(obj_meta$CNV)
obj_meta$CNV[obj_meta$CNV %in% c(2,5,6,7,8)] <- "C1"
obj_meta$CNV[obj_meta$CNV %in% c(3,4)] <- "C2"
obj_meta$CNV[obj_meta$CNV %in% c(1)] <- "C3"

## Prepare row annotation for complex heatmap
row_ha = rowAnnotation(CNV=obj_meta[,1], annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize=9), col=list(CNV=c("C1"="#1B9E77", "C2"="#D95F02", "C3"="#7570B3")), annotation_legend_param = list(CNV = list(at=names(table(obj_meta$CNV)), labels=paste0(names(table(obj_meta$CNV)), " (", c(82.4, 12.2, 5.4),"%", ")"))), show_legend=T)

## Plot Complex heatmap
pdf(paste0(sample, "_customized_infercnv.pdf"), width = 10.5, height = 6.5)
ht <- Heatmap(t(df), col=myColor, cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=F, show_row_dend=T, show_column_dend=F, show_heatmap_legend=T, left_annotation=row_ha, row_split=obj_meta$CNV, use_raster=T, name="Expression")
draw(ht, merge_legend = T)
for(i in 1:3) {
decorate_heatmap_body("Expression", {
    for ( i in end[1:22]) {
        grid.lines(c(i/end[23], i/end[23]), c(0, 1), gp = gpar(lty = 1, lwd = 1)) }}, row_slice = i)}
decorate_heatmap_body("Expression", {
    for ( j in 1:23) {
        grid.text(c(1:22, "X")[j], middle[j]/end[23], -0.02, rot=0, hjust=0.5, vjust=1, gp = gpar(fontsize=8))}}, row_slice = 3)
dev.off()

## Output cnv clone information for each cell
write.table(obj_meta[,1,drop=F], paste0(sample, "_cluster.txt"), row.names=T, col.names=F, sep="\t", quote=F)

