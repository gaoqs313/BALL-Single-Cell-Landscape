library(ggplot2)
library(dplyr)

## Load data (Supplementary Table 4)
df <- read.table("cnv_clone.txt", header = T, sep = "\t")

## Split by clone
dfs <- vector()
for(i in 1:nrow(df))
{
	if(grepl(";", df[i, 4]))
	{
		tmp <- as.character(unlist(strsplit(df[i,4], ";")))
		for(j in tmp)
		{
			dfs <- rbind(dfs, c(as.character(df[i, 1:3]), j))
		}
	}else
	{
		dfs <- rbind(dfs, as.character(df[i, ]))
	}
}
dfs <- as.data.frame(dfs)
colnames(dfs) <- colnames(df)
dfs$Blast <- as.integer(dfs$Blast)
dfs$Subclone <- as.integer(dfs$Subclone)

## Calculate Percetage of cells for subclons
dfs$Percentage <- round(dfs$Subclone/dfs$Blast*1000)/10

## Add clone number based on close size
dfs <- dfs[order(dfs$Sample, -dfs$Percentage), ]
nclone <- as.numeric(table(dfs$Sample))
tmp <- vector()
for(i in nclone)
{
        tmp <- c(tmp, 1:i)
}
dfs$Clone <- paste0("Clone", tmp)

## Sort the samples based on subclone size
dfs_sub1 <- as.data.frame(dfs %>% filter(Clone != "Clone1") %>% group_by(Sample) %>% summarise(Fraction = sum(Percentage)))
dfs_sub1 <- dfs_sub1[order(-dfs_sub1[,2]), ]
dfs_sub2 <- dfs[dfs$Percentage == 100, ]
dfs$Sample <- factor(dfs$Sample, levels = c(dfs_sub1[, 1], dfs_sub2[,1]), ordered = T)

## Specify clone color
cloneCol=c()
{
  cloneCol["Clone1"]='#1B9E77'
  cloneCol["Clone2"]='#D95F02'
  cloneCol["Clone3"]='#7570B3'
  cloneCol["Clone4"]='#E7298A'
  cloneCol["Clone5"]='#E6AB02'
  cloneCol["Clone6"]="#A6761D"
  cloneCol["Clone7"]="#666666"
  cloneCol["Clone8"]="#66A61E"
}

## Plot
pdf("Barplot_by_CNV_Subclone.pdf", width = 9, height = 5)
p <- ggplot(dfs, aes(x = Sample, y = Percentage, fill = Clone)) +
  geom_bar(stat="identity") +
  theme_bw() +
  facet_grid(~Subtype, scale='free_x', space="free") +
  theme(panel.spacing = unit(0.2, "lines")) +
  scale_y_continuous(limits = c(0,101), expand = c(0, 0)) +
  scale_fill_manual(values=cloneCol) +
  xlab("") +
  ylab("Percentage of subclones based on CNV") +
  theme(axis.title=element_text(size=12)) +
  theme(text = element_text(size=10, family = "Helvetica")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())
print(p)
dev.off()

