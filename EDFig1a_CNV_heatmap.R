library(ComplexHeatmap)
library(circlize)

## load input
input <- read.table("CopyNumber77S_4heatmap.tsv", header = T)

## convert to data frame by patients
inputs <- split(input, input[,1])
datain <- do.call(cbind, inputs)

## get ratio
ratio <- datain[, seq(7, ncol(datain), 8)]
rownames(ratio) <- datain[,8]
colnames(ratio) <- gsub(".mean", "", colnames(ratio))

## blacklist
bl <- c(2493:2494,3350:3420,4926:4928,7397:7407,8808:8818,12723:12724,13758:13766,22228:22240,23069:23080)
ratio[bl,] <- 0

## load subtype info
subtype <- read.table("sample_subtype.txt")
subtype <- subtype[order(subtype[,2], subtype[,1]), ]

## order
ratio <- ratio[, match(subtype[,1], colnames(ratio))]

## get telomere position
info <- as.data.frame(matrix(unlist(strsplit(rownames(ratio), ":|-")), ncol=3, byrow=T))
idxt <- c(1);
for(i in 2:nrow(info))
{
	if(info[i,1] != info[i-1,1])
	{
		idxt <- c(idxt, i)
	}
}
idxt <- c(idxt, nrow(info))
idxt <- idxt[2:24]

## get centromere location
info$V2 <- as.numeric(as.character(info$V2))
info$V3 <- as.numeric(as.character(info$V3))
centro <- read.table("hg19_centro.txt")
centro$V2 <- as.numeric(centro$V2)
idxc <- vector()
for(j in 1:nrow(centro))
{
        idxc <- c(idxc, which(as.character(info[,1])==as.character(centro[j,1]) & info[,2]<centro[j,2] & info[,3]>=centro[j,2]))
}

## get label position
label0 <- c(0, idxt[1:22])
label1 <- idxt[1:23]-label0
labelm <- label0 + label1/2
idxt <- idxt[1:22]

## color
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
  subtypeCol["Ph"]="magenta3"
  subtypeCol["Ph-like"]="brown"
  subtypeCol["PAX5alt"]="#FFA620"
  subtypeCol["TCF3-PBX1"]="darkgoldenrod4"
  subtypeCol["ZNF384"]="#A8DD00"
}

## prepare row annotation
row_ha = rowAnnotation(Subtype=subtype[,2], annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize=10), col=list(Subtype=subtypeCol))

## generate heatmap figure
myColor <- colorRamp2(seq(-1,1,0.02), colorRampPalette(colors = c("blue","white","red") )(101))
myChr <- c(1:22, 'X')

## Remove SJALL040066_D1
ratio <- ratio[, !colnames(ratio) %in% 'SJALL040066_D1']
row_ha = rowAnnotation(Subtype=subtype[subtype[,1] != 'SJALL040066_D1',2], annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize=10), col=list(Subtype=subtypeCol))

pdf("cnv_frequency_heatmap_v2.pdf", width=14, height=8)
ht_list <- Heatmap(as.matrix(t(ratio)), cluster_rows=F, cluster_columns=F, show_row_names=T, row_names_side="left", show_column_names=F, show_row_dend=F, show_column_dend=F, col=myColor, row_names_gp = gpar(fontsize=6), row_title=NULL, show_heatmap_legend=T, border=T, name="Log2(ratio)", heatmap_legend_param = list(direction = "vertical"), left_annotation=row_ha)
draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend=T, padding = unit(c(8, 2, 2, 2), "mm"))
decorate_heatmap_body("Log2(ratio)", {
    for ( i in idxt) {
        grid.lines(c(i/30375, i/30375), c(0, 1), gp = gpar(lty = 1, lwd = 1)) }
    for ( i in idxc) {
        grid.lines(c(i/30375, i/30375), c(0, 1), gp = gpar(lty = 2, col="grey60", lwd = 1)) }
    for ( j in 1:23) {
        grid.text(myChr[j], labelm[j]/30375, -0.012, rot=90, hjust=1, vjust=0.5) }
})
dev.off()

