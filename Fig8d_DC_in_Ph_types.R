
library(ggplot2)
library(dplyr)
library(reshape2)
library(Seurat)

## Load data
data <- read.table("../meta_for_all.txt", sep = "\t", header = T, row.names = 1)

## Add patient id to meta
data$Sample <- substr(rownames(data), 1, nchar(rownames(data))-19)

## Load subtype information
info <- read.table('/home/qgao/zdrive/Projects/scRNA_w_Ilaria/meta_20230105.txt', sep='\t', header=T)

## Add subtype to meta
data$Subtype <- info[match(data$Sample, info[, 3]), 8]

## Subset to Ph
df <- data[data$Subtype == 'Ph' & data$CellType != "Blast", ]

## Calculate DC percentage
perc <- vector()
for(i in unique(df$Sample))
{
	perc <- c(perc, sum(df$Sample == i & df$CellType == "DC") / sum(df$Sample == i))
}
ph <- data.frame(Sample = unique(df$Sample), Perc = perc*100)

## Add subtype
phtype <- read.table("../Ph_sub.txt")
ph$Type <- phtype[match(ph$Sample, phtype[,2]),1]
ph <- ph[order(ph$Perc, ph$Type), ]
ph$Sample <- factor(ph$Sample, levels = ph$Sample, ordered = T)

pdf("Percentage_by_DC_for_Ph.pdf", width=6, height=5)
p <- ggplot(ph, aes(x=Sample, y=Perc, fill=Type)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_y_continuous(limits = c(0,60), expand = c(0, 0)) +
  xlab("") +
  ylab("Percentage of DC out of all non-malignant cells") +
  theme(text = element_text(size=12), axis.title=element_text(size=12)) +
               theme(text = element_text(size=8, family = "Helvetica")) +
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = c(0.1, 0.8))
print(p)
dev.off()

###################################################################################################################

## Subset to DC
dc <- df[df$CellType == "DC", ]
dcg <- dc[dc$Andy %in% c('cDC1', 'cDC2', 'Pre-cDC', 'Pre-pDC', 'pDC'), ]

## Calculate pDC percentage
pperc <- vector()
ccount <- vector()
ttotal <- vector()
for(i in unique(dcg$Sample))
{
	ccount <- c(ccount, sum(dcg$Sample == i & dcg$Andy %in% c("pDC", "Pre-pDC")))
	ttotal <- c(ttotal, sum(dcg$Sample == i))
        pperc <- c(pperc, sum(dcg$Sample == i & dcg$Andy %in% c("pDC", "Pre-pDC")) / sum(dcg$Sample == i))
}
pph <- data.frame(Sample = unique(dcg$Sample), Total = ttotal, Number = ccount, Perc = pperc*100)
#          Sample Total    pDC      Perc
# SJBALL031128_D1   219    206  94.06393
#   SJPHALL007_D1   323    317  98.14241
#   SJPHALL010_D1    29     10  34.48276
#   SJPHALL020_D1    92     92 100.00000
#   SJPHALL020_R1    36     36 100.00000
#   SJPHALL021_D1   134    122  91.04478

###################################################################################################################

## Load data
merged <- readRDS("../Aggregate/Final_BALL_89S_Merged_Normal.rds")

## Subset to Ph
merged_ph <- subset(merged, subset = Subtype == "Ph")
merged_ph$Type <- phtype[match(merged_ph$orig.ident, phtype[,2]),1]
merged_ph$Type <- factor(merged_ph$Type, levels=c('Ph1', 'Ph2', 'Ph3'), ordered=T)
merged_ph$orig.ident <- factor(merged_ph$orig.ident, levels=ph$Sample, ordered = T)

## Select
selected <- c("CST3", 'IL3RA')
pdf("Ph_DC_markers.pdf", height=4.3, width=6)
p <- VlnPlot(merged_ph, features = selected, pt.size=0, group.by="orig.ident", split.by="Type", assay = "RNA", stack = T, flip = T, cols = c("#F8766D", "#00BA38", "#619CFF")) + NoLegend() + xlab("")
p <- p + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text = element_text(size=7, family = "Helvetica"))
print(p)
dev.off()

