library(Seurat)
library(dplyr)
library(magrittr)
library(data.table)
library(Rtsne)
library(readr)
library(ggplot2)
library(ggrepel)
library(sva)
library(uwot)
root="BALL-subtype-hg38/"
etc_root=paste0(root, "etc/")
source(paste0(etc_root, "init.R"))

## Load single cell pseudo bulk data
pseudo_rlog <- read.table("BALL_pseudobulk_rlog.tsv", header=T, row.names=1, sep='\t')

## Load reference bulk RNA-seq samples (GRCh38)
sampleInforStd_all=as.data.frame(fread(paste0(etc_root, "/input.STD4tSNE2046S.txt")))
sampleInforStd <- sampleInforStd_all
sampleNum=nrow(sampleInforStd)
for(i in 1:sampleNum){
  p=sampleInforStd$path[i]
  cat("Reading rds file: ", p, "\t", i, "/", sampleNum, "\n")
  expr_i=readRDS(paste0(p,".rds"))
  colnames(expr_i)=sampleInforStd$sample[i]
  if(i==1){
    rlogDF=expr_i
    geneIDs=rownames(rlogDF)
  }else{
    rlogDF=bind_cols(rlogDF, expr_i)
  }
}
rownames(rlogDF)=geneIDs

## Merge counts
rlogDF_subset <- rlogDF[match(rownames(pseudo_rlog), rownames(rlogDF)), ]
rlog_merged <- cbind(rlogDF_subset, pseudo_rlog)

## Merge meta
sampleInforPseudo <- data.frame(sample = colnames(pseudo_rlog), library = "SCRNA", subtype = "Others") 
sampleInfor <- rbind(sampleInforStd[, c(1:2,4)], sampleInforPseudo)

## Correct batch effect
modcombat = model.matrix(~1, data=sampleInfor)
rlog_bec = ComBat(dat=as.matrix(rlog_merged), batch=sampleInfor$library, mod=modcombat, par.prior=TRUE) %>% as.data.frame() %>% round(., digits = 2)

## Subset to pre-defined MAD genes
top1kGeneFile=paste0(etc_root, "STD-BALL-topMadGenes.txt");
top1kGenes=read_lines(top1kGeneFile)
rlog_forPlot <- rlog_bec[rownames(rlog_bec) %in% top1kGenes, ]

## Specify subtype colors
subtypeCol=c()
{
  subtypeCol["BCL2/MYC"]="seagreen2"
  subtypeCol["CDX2UBTF"]="#E52B50"
  subtypeCol["CRLF2_nonPhlike"]="#FFC0CB"
  subtypeCol["DUX4"]='grey40'
  subtypeCol["ETV6::RUNX1"]="gold2"
  subtypeCol["ETV6::RUNX1-like"]="deeppink"
  subtypeCol["HLF"]= "skyblue"
  subtypeCol["Hyperdiploid"]="#3E9F32"
  subtypeCol["iAMP21"]="lightslateblue"
  subtypeCol["IKZF1 N159Y"]="#E8E846"
  subtypeCol["KMT2A"]="#1F78B5"
  subtypeCol["KMT2A-like"]="grey60"
  subtypeCol["Low-hypodiploid"]="#1E90FF"
  subtypeCol["MEF2D"]="#66C2A6"
  subtypeCol["Near-haploid"]='#B8B828'
  subtypeCol["NUTM1"]='black'
  subtypeCol["Other"]='grey80'
  subtypeCol["PAX5 P80R"]="orangered"
  subtypeCol["PAX5alt"]='#FFA620'
  subtypeCol["BCR::ABL1"]="magenta3"
  subtypeCol["BCR::ABL1-like"]="brown"
  subtypeCol["TCF3::PBX1"]="darkgoldenrod4"
  subtypeCol["ZNF384"]="#A8DD00"
  subtypeCol["ZNF384-like"]="tomato3"
}

## Set parameters
N=1000 # top 1k gene
perp=30 # perplexity of tsne

## Plot for all pseudo-bulk samples
topMadDataTv=t(rlog_forPlot) %>% set_rownames(c())
set.seed(123)
rlogTsneOut = Rtsne(topMadDataTv, dims = 2, perplexity = perp, theta = 0.1, max_iter = 5000,
                        check_duplicates = F, partial_pca=T, num_threads = 4) # Run TSNE
tsneDF=sampleInfor %>% mutate(X=rlogTsneOut$Y[,1], Y=rlogTsneOut$Y[,2])
tsneDF$dataType <- "Bulk RNA"
tsneDF$dataType[tsneDF$subtype == "Others"] <- "scRNA pseudobulk"

tsneGP=ggplot() + xlab("tSNE-1") + ylab("tSNE-2")+
  geom_point(data=subset(tsneDF, dataType  == "Bulk RNA"), aes(X, Y, color=subtype, shape=dataType), size=0.2)+
  geom_point(data=subset(tsneDF, dataType  == "scRNA pseudobulk" ), aes(X, Y, color=subtype, shape=dataType), size=1, stroke = 0.8)+
  theme_bw() +
  scale_shape_manual(values = c(16, 1), labels = c("Bulk RNA (2046)", "scRNA pseudobulk (89)"), name = "Data type") +
  scale_color_manual(values = subtypeCol, name = "BALL subtype") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(text = element_text(family = "Helvetica"),
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(tsneGP, width=8.5, height=5, filename="tSNE_pseudobulk.pdf", useDingbats=FALSE)


