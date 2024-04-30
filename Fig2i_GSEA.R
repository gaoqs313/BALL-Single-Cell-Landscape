library(ggplot2)
library(viridis)

## Load GSEA result
input <- read.table("GSEA_manual.txt", header=T)
input$LFDR <- input$FDR
input$LFDR[input$LFDR==0] <- 1e-6
input$LFDR <- -log10(input$LFDR)
input$Pathway <- gsub("HALLMARK_", "", input$Pathway)

## Bubble plot
pdf("Pathway_bubble.pdf", width = 6, height = 5.6)
p <- ggplot(data=input, aes(x=Subtype, y=Pathway)) +
  geom_point(aes(col = NES, size = LFDR)) +
  scale_color_continuous(name = "NES", type = "viridis",  breaks = c(1, 2, 3, 4)) +
  scale_size_continuous(name = "FDR", breaks = rev(c(-log10(c(0.01, 0.0001, 0.000001)))), labels = rev(c("1e-2", "1e-4", "1e-6"))) +
  theme_bw() +
  facet_grid(Program~., scales="free", space="free") +
  theme(strip.background=element_rect(fill="white")) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=8), axis.text.y = element_text(size=8)) +
  theme(legend.position = "left")
print(p)
dev.off()

