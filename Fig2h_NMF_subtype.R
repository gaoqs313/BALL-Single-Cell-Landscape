library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer)

subtype <- "DUX4"

## Specify parameters
myDir <- 'cNMF_output'
cor_min <- 0
cor_max <- 0.8

## Load data
kdt <- read.table(paste0(subtype, "_k_dt.txt"), header = F)
all.score.df <- data.frame()
for (i in 1:nrow(kdt))
{
        score.df <- read.table(paste0(myDir, "/", kdt$V1[i],"/program.Zscore.txt"), header = T,sep = "\t",stringsAsFactors = F, check.names=F)
        if (i==1) { all.score.df=score.df }
        if (i!=1) { all.score.df=all.score.df%>%inner_join(score.df,by="gene") }
}

## Clean up
rownames(all.score.df) <- all.score.df$gene
all.score.df$gene <- NULL
all.score.df <- all.score.df[rowSums(is.na(all.score.df)) == 0,]

## Remove bacground programs
if(file.exists(paste0(myDir, "/", subtype, "/background_program.txt")))
{
    maybe.bg <- read.table(paste0(myDir, "/", subtype, "/background_program.txt"))
    maybe.bg <- as.character(maybe.bg[,1])
    all.score.rm.df <- all.score.df[,setdiff(colnames(all.score.df),maybe.bg)]
}else
{
    all.score.rm.df <- all.score.df
}

## Pearson correlation
all.score.rm.df.cor <- cor(all.score.rm.df,method = "pearson")
all.score.rm.df.cor[all.score.rm.df.cor < cor_min]=cor_min
all.score.rm.df.cor[all.score.rm.df.cor > cor_max]=cor_max

## Prepare heatmap color
myColor <- colorRamp2(seq(0, 0.8, 0.008), colorRampPalette(colors = c("white", "yellow", "red", "#67001F"))(101))

## Plot heatmap
pdf(paste0(subtype, "_NMF_correlation_heatmap.pdf"), width=8)
ht <- Heatmap(as.matrix(all.score.rm.df.cor), cluster_columns=T, col=myColor, cluster_rows=T, show_row_dend=F, show_column_dend=F, show_row_names=F, show_column_names=F, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6), name="Correlation", show_heatmap_legend = T)
draw(ht, merge_legend=T)
decorate_heatmap_body("Correlation", {
    grid.lines(c(0, 0), c(51/61, 1), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(10/61, 10/61), c(51/61, 1), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(0, 10/61), c(51/61, 51/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(0, 10/61), c(1, 1), gp = gpar(lty = 1, lwd = 3, col = 1))

    grid.lines(c(10/61, 10/61), c(41/61, 51/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(20/61, 20/61), c(41/61, 51/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(10/61, 20/61), c(41/61, 41/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(10/61, 20/61), c(51/61, 51/61), gp = gpar(lty = 1, lwd = 3, col = 1))

    grid.lines(c(20/61, 20/61), c(31/61, 41/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(30/61, 30/61), c(31/61, 41/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(20/61, 30/61), c(31/61, 31/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(20/61, 30/61), c(41/61, 41/61), gp = gpar(lty = 1, lwd = 3, col = 1))

    grid.lines(c(30/61, 30/61), c(21/61, 31/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(40/61, 40/61), c(21/61, 31/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(30/61, 40/61), c(21/61, 21/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(30/61, 40/61), c(31/61, 31/61), gp = gpar(lty = 1, lwd = 3, col = 1))

    grid.lines(c(40/61, 40/61), c(0/61, 21/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(61/61, 61/61), c(0/61, 21/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(40/61, 61/61), c(0/61, 0/61), gp = gpar(lty = 1, lwd = 3, col = 1))
    grid.lines(c(40/61, 61/61), c(21/61, 21/61), gp = gpar(lty = 1, lwd = 3, col = 1))
})
dev.off()

