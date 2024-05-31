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

==========================================================
## Establish a normal human B-cell development atlas

We provide a notebook outlining our parameters to [construct the cross-ontogeny B-cell development atlas](https://htmlpreview.github.io/?https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig3.0_BDev_ProjectionSetup.nb.html) from the component datasets and a separate notebook to [characterize human B cell development at the single cell level](https://htmlpreview.github.io/?https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig3.1_BDevelopment_Atlas_Characterization.nb.html). 

Notebooks and scripts related to analyses in Figure 3 
```
Fig3.0_BDev_ProjectionSetup.Rmd
Fig3.1_BDevelopment_Atlas_Characterization.Rmd
Fig3f_CEBPA_ATAC.R
```

==========================================================
## B-ALL projection along B cell development

We provide a notebook outlining how we performed [scRNA-seq projections of B-ALL transcriptomes](https://htmlpreview.github.io/?https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig4.0_runBALL_projections_updated.nb.html) onto normal B cell development and use these projection results to perform [B-ALL composition analysis](https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig4.1_BALL_composition_Analysis.ipynb). We summarise the composition analyses results in a complex heatmap depicting [B-ALL compositional clusters](https://htmlpreview.github.io/?https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig4.2_BALL_CompositionHeatmap.nb.html).

Notebooks and scripts related to analyses in Figure 3 
```
Fig4.0_runBALL_projections_updated.Rmd
Fig4.1_BALL_composition_Analysis.ipynb
Fig4.2_BALL_CompositionHeatmap.Rmd
```

==========================================================
## Bulk RNA-seq analysis using LASSO regression

Bulk RNA-seq analysis using LASSO regression (ED Fig. 7)
```
ED Fig. 7_DevState_LASSOregression.Rmd
ED Fig. 7_DevState_LASSOregression.html
```

B-ALL Dev State Analysis (Fig. 5a-c, ED Fig. 6-8)
```
Fig5a-c_BALL_Multipotency_Characterization.ipynb
```

==========================================================
## B-ALL Composition Multipotency Score

B-ALL Composition Multipotency Score (Fig. 5 d-l, Fig. 6, ED Fig. 9-10)
```
Fig5i-l&6_BALL_MultipotencyScore.Rmd
Fig5i-l&6_BALL_MultipotencyScore.html
```

