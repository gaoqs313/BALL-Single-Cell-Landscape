
library(ggplot2)
library(dplyr)
library(Seurat)

## Load data
df <- read.table("Cell_count_in_different_types.txt", sep = "\t", header = T, row.names = 1)

## Merge Ph-like
df[df == "Ph-like_CRLF2"] <- "Ph-like"
df[df == "Ph-like_non CRLF2"] <- "Ph-like"

## Keep Monocyte >= 10
mono <- df[df$Monocyte >= 10, ]
#57

## Calculate CD16 percentage
ratio <- data.frame(cd16 = mono$CD16.Mono/mono$Monocyte*100, Subtype = mono$Subtype)
ratio$Sample <- rownames(mono)

## Output
write.table(ratio, "SuppTable_for_Fig8b.txt", quote=F, sep="\t", row.names=F, col.names=T)

## Specify subtype colors
subtypeCol=c()
{
  subtypeCol["DUX4"]='grey40'
#  subtypeCol["ETV6-RUNX1"]="gold2"
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

## Calculate median
ratio_median <- as.data.frame(ratio %>% group_by(Subtype) %>% summarise(cd16 = median(cd16)))
ratio$Subtype <- factor(ratio$Subtype, levels = ratio_median[order(ratio_median[,2]), 1], ordered = T)

pdf("Percentage_CD16.pdf", width = 8, height = 6)
p <- ggplot(ratio, aes(x=Subtype, y=cd16, color=Subtype)) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(position = position_jitter(width=0.2, height=0),size=3) +
	scale_color_manual(values=subtypeCol) +
	xlab("") + ylab("Percentage of CD16+ cells in monocytes") +
	theme_bw() +
	theme(legend.position = "none") +
	theme(text = element_text(size=14, family = "Helvetica")) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
print(p)
dev.off()


