library(ComplexHeatmap)

################################################
##### Oncoprint using mutation & meta data #####
################################################

mm_input <- read.table("Input_mutation_and_meta.tsv", sep = "\t", header = T, row.names = 1, check.names = F)
driver <- read.table("Driver.tsv", sep = "\t")
driver <- driver[match(rownames(mm_input)[1:46], driver[,1]), ]
driver[driver == "Cell cycle/tumor suppressor"] <- "Cell cycle"
driver[driver == 'DNA damage'] <- "Others"
driver[driver == 'Drug response'] <- "Others"
driver[,2] <- factor(driver[,2], levels=c("B cell", "Cell cycle", "Signaling", "Transcription", "Epigenetic", "Others"), ordered=T)

##################
## Prepare data ##
##################

## extract data
mut <- mm_input[1:46, ]
meta <- as.data.frame(t(mm_input[c(47:51), ]))
meta$WGS <- factor(meta$WGS, levels=c('YES', 'NO'), ordered=T)
meta$WES <- factor(meta$WES, levels=c('YES', 'NO'), ordered=T)
meta$RNASeq <- factor(meta$RNASeq, levels=c('YES', 'NO'), ordered=T)

################################
## Prepare color and function ##
################################

## define colors for mutations
colMUT = c(SNV = "blue2",
           INDEL = "goldenrod2", 
           AMP = "#E41A1C",
           DEL = "#4DAF4A")

## define colors for meta data
colors = list("Subtype" = c("BCR::ABL1" = "magenta3", 
                             "BCR::ABL1-like" ="brown", 
                             "DUX4-R" = "grey40",
                             "ETV6::RUNX1-like" = "deeppink",
			     "Hyperdiploid" = "#3E9F32",
			     "iAMP21" = "lightslateblue",
			     "KMT2A-R" = "#1F78B5",
			     "Low-hypodiploid" = "#1E90FF",
			     "MEF2D-R" = "#66C2A6",
			     "Near-haploid" = '#B8B828',
			     "Other" = 'grey75',
			     "PAX5alt" = "#FFA620",
			     "TCF3::PBX1" = "darkgoldenrod4",
			     "ZNF384-R" = "#A8DD00"),
              "WGS" = c("YES" = "#1F77B4",
                        "NO" = "#B8B8B8"),
              "WES" = c("YES" = "#1F77B4",
                        "NO" = "#B8B8B8"),
              "RNASeq" = c("YES" = "#1F77B4",
                            "NO" = "#B8B8B8"),
              "Stage" = c("Initial diagnosis" = "#1F77B4",
                            "Recurrence" = "#B8B8B8"))

## define alternation function
alter_fun = function(x, y, w, h, v) {
  n = sum(v)
  h = h*0.9
  grid.rect(x, y, w-unit(0.2,'mm'), h-unit(0.2, 'mm'), gp=gpar(fill="grey90", col=NA))
  if(n) grid.rect(x, y-h*0.5+1:n/n*h, w-unit(0.2, 'mm'), 1/n*h, gp=gpar(fill=colMUT[names(which(v))], col=NA), just='top')
}

#########################
## Prepare annotations ##
#########################

top_annotation = HeatmapAnnotation(Subtype = meta$Subtype,
			    Stage = meta$Stage,
			    WGS = meta$WGS,
                            WES = meta$WES,
                            RNASeq = meta$RNASeq,
                            which = c("column"),
                            simple_anno_size = unit(2, "mm"),
                            col = colors,
                            border= TRUE,
                            na_col = '#B8B8B8', 
                            annotation_name_gp = gpar(fontsize = 7),
                            annotation_name_side = "left")

####################
## Plot oncoprint ##
####################

## change row names
cnt <- apply(!is.na(mut), 1, sum) 
mutn <- mut
rownames(mut) <- paste0(rownames(mut), " (", cnt, ")")

oncoprint = oncoPrint(mut, get_type = function(x) strsplit(x, ";")[[1]],
                alter_fun = alter_fun, col = colMUT, 
                row_names_gp = gpar(fontsize = 6), 
		column_title_gp = gpar(fontsize = 8),
		row_title_side = "right", 
		row_title_gp = gpar(fontsize = 8),
		column_title_rot = 90,
                show_column_names = F,
                show_pct = F, 
                row_names_side = "left",
                row_order = rownames(mut),
                right_annotation = NULL, 
                column_split = meta$Subtype,
	        row_split = driver[,2],
		row_gap = unit(0.5, "mm"),
                top_annotation = top_annotation)

pdf("Oncoprint_scRNA.pdf", width=7.2, height=8.8)
draw(oncoprint, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend=T, padding = unit(c(6, 2, 10, 2), "mm"))
dev.off()

