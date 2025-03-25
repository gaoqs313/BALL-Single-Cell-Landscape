library(reshape2)
library(ggplot2)

## Load NMF signature genes from each subtype
subtype <- c('DUX4', 'ETV6-RUNX1-like', 'Hyperdiploid', 'KMT2A', 'MEF2D', 'Near-haploid', 'PAX5alt', 'Ph', 'Ph-like', 'TCF3-PBX1', 'ZNF384')
results <- vector()
for(i in 1:5)
{
	topGene <- vector()
	## Loop through subtypes
	for(s in subtype)
	{
		gtop <- read.table(paste0("../A.Subtype/", s, "/", s, "_top30_genes_for_meta_heatmap.txt"))
		topGene <- c(topGene, gtop[(30*(i-1)+1):(30*i), 1])
	}
	counts <- table(topGene)
	topGene_replaced <- as.data.frame(table(as.vector(counts[topGene])))
	if(i == 1)
	{
		results = topGene_replaced
	}else
	{
		results = merge(results, topGene_replaced, by = "Var1", all = T)
	}
}
results[is.na(results)] <- 0
colnames(results) <- c('N', 'A', 'B', 'C', 'D', 'E')
df <- melt(results)

## Barplot
pdf("Signature_specificity.pdf", height=5, width = 8)
p <- ggplot(df, aes(x=N, y=value, fill=N)) +
        geom_bar(stat = "identity", colour="black") +
        theme_bw() +
	facet_wrap(.~variable, ncol = 5) +
        xlab("Number of subtypes") +
        ylab("Number of genes") +
        theme(axis.text.x = element_text(hjust = 0.5, size = 8, family = "Helvetica")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(text = element_text(size=14, family = "Helvetica")) +
	theme(legend.position = "none")
print(p)
dev.off()
