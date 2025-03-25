
library(Seurat)
library(ggplot2)
library(reshape2)

## Load data
df_meta <- read.table("/home/qgao/zdrive/Projects/scRNA_w_Ilaria/Aggregate/Final_BALL_89S_Merged_Metadata.txt", header = T, sep = "\t", row.names = 1)

## Remove Other
df_meta <- df_meta[df_meta$CellType != "Blast", ]

## Specify cell type color
celltypeCol=c()
{
  celltypeCol["B"]='#4DAF4A'    #green
  celltypeCol["Blast"]='#377EB8' #blue
  celltypeCol["DC"]='#FF7F00'    #orange
  celltypeCol["Erythroid"]='#E41A1C'  #red
  celltypeCol["HSC/MPP and pro."]='#00CDCD' #dark cyan
  celltypeCol["Monocyte"]='#984EA3' #purple
  celltypeCol["Plasmablast"]='#EBBC2E' #gold
  celltypeCol["T_NK"]='#A65628'    #brown
}

## Get count
count <- as.matrix(table(df_meta[,c(1, 7)]))

## Calculate percentage
perc <- count / apply(count, 1, sum) * 100

## Output
write.table(perc, "Supplementary_table.txt", row.names=T, col.names=NA, sep = "\t", quote=F)

## Prepare data
df <- melt(perc)
meta <- unique(df_meta[, c(1,5)])
df$Subtype <- meta[match(df[,1], meta[,1]), 2]

pdf("Percentage_by_normal.pdf", width=9, height=5)
p <- ggplot(df, aes(x=orig.ident, y=value, fill=CellType)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_grid(~Subtype, scale='free_x', space="free") +
  theme(panel.spacing = unit(0.2, "lines")) +
  scale_y_continuous(limits = c(0,101), expand = c(0, 0)) +
  scale_fill_manual(values = celltypeCol) + 
  xlab("") +
  ylab("Percentage of non-malignant cell types") +
  theme(legend.position = "right") +
  theme(axis.title=element_text(size=12)) +
  theme(text = element_text(size=10, family = "Helvetica")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(p)
dev.off()
