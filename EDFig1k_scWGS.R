
library(zoo)
library(ggplot2)

## Load samples
sample <- read.table("id.txt")

## Define parameter
myDir <- "/home/qgao/zdrive/Projects/scRNA_w_Ilaria/TCF3-PBX1/Revision/hg38/SJE2A067"
#myNorm <- c("A11", "B11", "C11", "D11", "E11", "F11", "G11", "H11")
myNorm <- c("B11", "C11", "E11", "F11", "G11", "H11")

## Load genome
genomeSize <- read.table("hg38_length.txt", comment.char = "", row.names = 1)
colnames(genomeSize) <- c("chrSize", "Color")
genomeSize$prevLen <- cumsum(c(0,genomeSize$chrSize[1:21]+6e6))
genomeSize$plotPos <- genomeSize$prevLen + genomeSize$chrSize / 2

## Load data
myList <- list()
for(i in sample[,1])
{
	df <- read.table(paste0(myDir, "/", i, "/coverage.bin.segments.txt"), header = T)
	
	## Remove X Y
	df <- df[!df$rname %in% c('X', 'Y'), ]

	## Add plot position
	df$prevLen <- genomeSize$prevLen[match(df$rname, rownames(genomeSize))]
	df$plotSTART <- df$startpos + df$prevLen
	df$plotEND <- df$endpos + df$prevLen
	df$prevLen <- NULL

	## Normalize by global median
	df$selfnormDepth <- df$meandepth / median(df$meandepth, na.rm = T)
	
	## Add to list
	myList[[i]] <- df
}

## Generate normalization factor using normal cells
df_norm <- vector()
for(i in myNorm)
{
	df_norm <- cbind(df_norm, myList[[i]]$selfnormDepth)
}
normFactor <- apply(df_norm, 1, median)

## Normlize using normal cells
for(i in sample[,1])
{
	myList[[i]]$meanRatio <- normFactor
	myList[[i]]$NormDepth <- myList[[i]]$meandepth / normFactor
	myList[[i]]$CopyNumber <- 2 * myList[[i]]$NormDepth / median(myList[[i]]$NormDepth)
	## Output
	write.table(myList[[i]], paste0(myDir, "/", i, "/scWGS.CNVcalling.txt"), row.names=F, col.names=T, sep="\t", quote=F)
}

## Plot
pdf("Coverage_segment.pdf", width = 20, height = 3)
for(i in sample[,1])
{
	dfp <- myList[[i]]
	dfp <- dfp[complete.cases(dfp), ]
	dfp$rname <- as.factor(dfp$rname)
	p <- ggplot(dfp) +
	geom_segment(aes(x = plotSTART, y = CopyNumber, xend = plotEND, yend = CopyNumber, color = rname)) +
    	theme_bw() +
    	ylim(-1, round(max(dfp$CopyNumber)) + 2) +
        xlab("") + ylab("Copy number") +
    	scale_color_manual(values = genomeSize$Color) +
        scale_x_continuous(limits = c(0,3100000000), breaks = genomeSize$plotPos, labels = 1:22, expand = c(0, 0)) +
    	theme(legend.position = "none") +
	ggtitle(i) +
    	theme(text = element_text(family = "Helvetica"),
          axis.title=element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
	print(p)
}
dev.off()



