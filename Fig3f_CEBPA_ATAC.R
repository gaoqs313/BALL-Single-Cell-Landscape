library(EnrichedHeatmap)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)

samp <- read.table("BCell_ATAC.txt")

counts <- vector()

for(i in 1:nrow(samp))
{
  df <- fread(paste(samp[i,1], ".wig", sep=""))
  dfg <- GRanges(seqnames=df[[1]], ranges=IRanges(df[[2]]+1, df[[3]]), count=df[[4]])
  window=makeWindows(dfg, w = 50)
  mtch = as.matrix(findOverlaps(window, dfg))
  window$count <- dfg$count[mtch[,2]]
  window <- window[start(window)>=33258000 & end(window)<=33305000, ]
  win <- as.data.frame(window)
  counts <- cbind(counts, win$count)
}

rownames(counts) <- apply(win[, 1:3], 1, paste, collapse="_")
colnames(counts) <- samp[,1]

data <- as.data.frame(apply(counts, 2, function(x) smooth.spline(x)$y))
data$x <- 1:nrow(data)
data[data < 0] <- 0

myList <- list()

#LMPP
lmpp <- data[, c(1:3, 26)]
lmpp$Z <- apply(data[, 1:3], 1, mean)
lmpp <- melt(lmpp, id.vars="x")
p <- ggplot(lmpp, aes(x=x, y=value, fill=variable)) +
  geom_area(position="identity") +
  scale_fill_manual(values=c(rep(adjustcolor("#45B500", alpha.f=0.05), 3), adjustcolor("#45B500", alpha.f=0.2))) +
  geom_line(aes(color=variable)) +
  scale_color_manual(values=c(rep(adjustcolor("#45B500", alpha.f=0.2), 3),"#45B500")) +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(breaks=c(0, 3), limits=c(0, 3)) +
  annotate("text",  x=460, y=2.5, label = "LMPP (n=3)", colour = "#45B500", fontface="bold", size = 6) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(text = element_text(size=20, family = "Helvetica")) +
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0.2,0.2,-0.6,0.2), "cm")) +
  xlab("") + ylab("")
myList[[1]] <- p

#MLP
mlp <- data[, c(4:6, 26)]
mlp$Z <- apply(data[, 4:6], 1, mean)
mlp <- melt(mlp, id.vars="x")
p <- ggplot(mlp, aes(x=x, y=value, fill=variable)) + 
  geom_area(position="identity") + 
  scale_fill_manual(values=c(rep(adjustcolor("#00C087", alpha.f=0.05), 3), adjustcolor("#00C087", alpha.f=0.2))) +
  geom_line(aes(color=variable)) +
  scale_color_manual(values=c(rep(adjustcolor("#00C087", alpha.f=0.05), 3),"#00C087")) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(breaks=c(0, 3), limits=c(0, 3)) +
  annotate("text",  x=460, y=2.5, label = "MLP (n=3)", colour = "#00C087", fontface="bold", size = 6) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(text = element_text(size=20, family = "Helvetica")) +
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0.2,0.2,-0.6,0.2), "cm")) +
  xlab("") + ylab("")
myList[[2]] <- p

#CLP
clp <- data[, c(7:12, 26)]
clp$Z <- apply(data[, 7:12], 1, mean)
clp <- melt(clp, id.vars="x")
p <- ggplot(clp, aes(x=x, y=value, fill=variable)) + 
  geom_area(position="identity") + 
  scale_fill_manual(values=c(rep(adjustcolor("#F8766D", alpha.f=0.05), 6), adjustcolor("#F8766D", alpha.f=0.2))) +
  geom_line(aes(color=variable)) +
  scale_color_manual(values=c(rep(adjustcolor("#F8766D", alpha.f=0.05), 6),"#F8766D")) +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(breaks=c(0, 3), limits=c(0, 3)) +
  annotate("text",  x=460, y=2.5, label = "CLP (n=6)", colour = "#F8766D", fontface="bold", size = 6) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(text = element_text(size=20, family = "Helvetica")) +
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0.2,0.2,-0.6,0.2), "cm")) +
  xlab("") + ylab("")
myList[[3]] <- p

#Pre-Pro-B
preprob <- data[, c(13:16, 26)]
preprob$Z <- apply(data[, 13:16], 1, mean)
preprob <- melt(preprob, id.vars="x")
p <- ggplot(preprob, aes(x=x, y=value, fill=variable)) + 
  geom_area(position="identity") + 
  scale_fill_manual(values=c(rep(adjustcolor("#9C8DFF", alpha.f=0.05), 4), adjustcolor("#9C8DFF", alpha.f=0.2))) +
  geom_line(aes(color=variable)) +
  scale_color_manual(values=c(rep(adjustcolor("#9C8DFF", alpha.f=0.05), 4),"#9C8DFF")) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(breaks=c(0, 3), limits=c(0, 3)) +
  annotate("text",  x=460, y=2.5, label = "Pre-Pro-B (n=4)", colour = "#9C8DFF", fontface="bold", size = 6) +
  theme(plot.title = element_text(color="#9C8DFF", face="bold")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(text = element_text(size=20, family = "Helvetica")) +
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0.2,0.2,-0.6,0.2), "cm")) +
  xlab("") + ylab("                             ATAC-seq (rpm/bp)")
myList[[4]] <- p

#ProB
prob <- data[, c(17:22, 26)]
prob$Z <- apply(data[, 17:22], 1, mean)
prob <- melt(prob, id.vars="x")
p <- ggplot(prob, aes(x=x, y=value, fill=variable)) +
  geom_area(position="identity") +
  scale_fill_manual(values=c(rep(adjustcolor("#FF61C7", alpha.f=0.05), 6), adjustcolor("#FF61C7", alpha.f=0.2))) +
  geom_line(aes(color=variable)) +
  scale_color_manual(values=c(rep(adjustcolor("#FF61C7", alpha.f=0.05), 6),"#FF61C7")) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(breaks=c(0, 3), limits=c(0, 3)) +
  annotate("text",  x=460, y=2.5, label = "Pro-B (n=6)", colour = "#FF61C7", fontface="bold", size = 6) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(text = element_text(size=20, family = "Helvetica")) +
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0.2,0.2,-0.6,0.2), "cm")) +
  xlab("") + ylab("")
myList[[5]] <- p

pdf("CEBPA_Bdev.pdf", width = 6, height = 5)
ggarrange(plotlist=myList, ncol=1, nrow=5)
dev.off()
