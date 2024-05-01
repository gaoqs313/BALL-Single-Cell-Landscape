# BALL-Single-Cell-Landscape

This is the repository for scripts used in our manuscript "SINGLE CELL DISSECTION OF DEVELOPMENTAL ORIGINS AND TRANSCRIPTIONAL HETEROGENEITY IN B-CELL ACUTE LYMPHOBLASTIC LEUKEMIA"


==========================================================
## Merge 89 SC samples and generate UMAP

Generate merged UMAP color-coded by subtype (Fig. 1b), cell type (ED Fig. 1c), or tissue (ED Fig. 1d)
```
Rscript Fig1b_UMAP_Merge.R
```

Generate tSNE of bulk RNA-seq and pseudo-bulk RNA-seq samples (ED Fig. 1b)
```
Rscript EDFig1b_tSNE_bulk.R
```

Check somatic mutations in this cohort (ED Fig. 1a)
```
Rscript EDFig1a_oncoprint.R
```

==========================================================
## Dissect copy number heterogeneity using inferCNV

Run inferCNV for each sample
```
Rscript Fig2_infercnv_command.R
```

Plot copy number clone distribution (Fig. 2a)
```
Rscript Fig2a_CNV_Barplot.R
```

Optimize inferCNV output using Complex Heatmap (Fig. 2b, Fig. 2d, ED Fig. 2c, ED Fig. 2e-g)
```
Rscript Fig2bd_CNV_custom.R
```

Plot UMAP for each sample color-coded by cell type, copy number clone, or fusion status (Fig. 2c, Fig. 2e, ED Fig. 2d and ED Fig. 2h)
```
Rscript Fig2ce_UMAP.R
```

Validate copy number clone using scWGS (Fig. 2f and ED Fig. 2i)
```
Rscript Fig2f_scWGS.R
```

==========================================================
## Dissect transcriptomic heterogeneity using cNMF

Run cNMF for each sample
```
Rscript Fig2_cnmf_command.sh
```

Plot pairwise correlation heatmap for each subtype (Fig. 2h and ED Fig. 3a)
```
Rscript Fig2h_NMF_subtype.R
```

Use GSEA to define each NMF program (Fig. 2i)
```
Fig2i_GSEA.R 
```

Check correlation of NMF programs across subtypes (Fig. 2j)
```
Rscript Fig2j_NMF_pairwise.R
```

Check subtype specificity of NMF signature genes (ED Fig.3b)
```
Rscript EDFig3b_NMF_specificity.R
```

