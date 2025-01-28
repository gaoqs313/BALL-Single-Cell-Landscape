
library(ggplot2)

## Format
module_score$Cell <- rownames(module_score)
df_score <- melt(module_score, id = c("orig.ident", "Subtype", "Cell"))

## Calculate median
median_score <- as.data.frame(module_score %>% group_by(orig.ident) %>% summarise(Escore = median(E)))
median_score <- median_score[order(median_score[,2]), ]

## Add DUX4 type
dtype <- read.table("DUX4_sub.txt", row.names=1)
df_score$Type <- dtype[match(df_score[,1], rownames(dtype)),3]
df_score$orig.ident <- factor(df_score$orig.ident, levels=median_score[,1], ordered=T)

pdf("Violin_program_score.pdf", width=6)
p <- ggplot(df_score, aes(x=orig.ident, y=value, fill=Type)) +
  geom_violin() +
  geom_boxplot(width=0.1, outlier.color=NA)+
  theme_bw() +
  facet_grid(variable~.) +
  xlab("") +
  ylab("Program score") +
  theme(legend.key = element_blank(), strip.background = element_rect(fill="white"))+
  theme(axis.text.x = element_text(size = 6, family = "Helvetica", angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top")
print(p)
dev.off()

