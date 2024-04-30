library(ggplot2)

## Load data (Supplementary Table 4)
data <- read.table("scWGS_segment_copy.txt", header = T)

## Sort samples by 1ql
tc <- data[data$Frag == '1ql', ]
tco <- tc[order(tc$CopyNumber), ]
data$Sample <- factor(data$Sample, levels = tco$Sample, ordered = T)

## Add gain/loss based on segment
data$Status <- "Normal"
data$Status[data$CopyNumber > 2.5] <- "Gain"
data$Status[data$CopyNumber < 1.5] <- "Loss"
data$Frag <- factor(data$Frag, levels = c("1p", "1ql", "1qr", "9p", "9q", "19p"), ordered=T)

## Create fake data to control y-axis limit for each region
fake1 <- data[data$Sample == "A01", ]
fake1$CopyNumber <- 0
fake2 <- data[data$Sample == "A01", ]
fake2$CopyNumber <- c(3,5.2,5.2,3,5,3)
fake <- rbind(fake1, fake2)
fake$Sample <- NA

## Plot
pdf("SJE2A063_Dotplot.pdf", width=14, height=8, useDingbats = F)
p <- ggplot(data=data, aes(x=Sample, y=CopyNumber)) +
	geom_point(aes(color = Status)) +
	geom_point(data=fake, x=NA) +
	xlab("Cells") +
	ylab("Copy number") +
	scale_color_manual(values=c("red", "blue", "black")) +
	facet_grid(Frag ~ Clone, scales = "free", space = "free_x") +
	geom_hline(yintercept = 2, col = "grey30", linetype = "dashed") +
	theme_bw() +
	theme(text = element_text(size=16, family = "Helvetica")) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	theme(strip.background =element_rect(fill="white"))
print(p)
dev.off()

