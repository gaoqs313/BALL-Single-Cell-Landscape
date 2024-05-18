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

Establish a single-cell atlas (Fig. 3c)
```
Fig3c_BDev_ProjectionSetup.Rmd
Fig3c_BDev_ProjectionSetup.html
```

==========================================================
## B-ALL projection along B cell development

Project B-ALL cells onto B cell development reference (Fig. 4a, Fig. 4c)
```
Fig4ac_runBALL_projections.Rmd
Fig4ac_runBALL_projections.html
```

B-ALL composition analysis (Fig. 4d)
```
Fig4d_BALL_CompositionAnalysis.Rmd
Fig4d_BALL_CompositionAnalysis.html
```

ADD SCRIPTS FOR FIGURE 4 E-H UMAPS

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

