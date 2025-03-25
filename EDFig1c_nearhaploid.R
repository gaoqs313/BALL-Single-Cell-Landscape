
library(ggplot2)

## Load data
df <- read.table("nearhaploid.txt", header = F, sep = "\t")
colnames(df) <- c("Sample", "Chr", "Percentage", "Status")

## Sort
df$Chr <- factor(df$Chr, levels = c(1:22, 'X'))
df$Sample <- gsub("_D1|_D", "", df$Sample)
df$Status <- factor(df$Status, levels = c("Partial", "Complete"), ordered = T)

## Plot
pdf("Barplot_nearhaploid.pdf", width = 14, height=6)
p <- ggplot(df, aes(x = Chr, y = Percentage, fill = Status)) +
  geom_bar(stat="identity") +
  facet_grid(. ~ Sample, scales = "free", space = "free") +
  theme_bw() +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  scale_fill_manual(values = c("#00b4d8", "darkblue")) +
  ylab("Percentage of cells with copy number loss") +
  xlab("") +
  theme(text = element_text(size=10, family = "Helvetica")) +
  theme(strip.background=element_rect(fill="white"), strip.text.x = element_text(size = 7)) +
  theme(axis.title=element_text(size=14)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()

