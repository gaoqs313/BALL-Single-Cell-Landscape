# BALL-Single-Cell-Landscape

This is the repository for scripts used in our manuscript "SINGLE CELL DISSECTION OF DEVELOPMENTAL ORIGINS AND TRANSCRIPTIONAL HETEROGENEITY IN B-CELL ACUTE LYMPHOBLASTIC LEUKEMIA"


==========================================================
## Merge 89 SC samples and generate UMAP

Generate merged UMAP color-coded by subtype (Fig. 1b), cell type (ED Fig. 1c), or tissue (ED Fig. 1d)
```
Rscript Fig1c-d_UMAP_Merge.R
```

Generate tSNE of bulk RNA-seq and pseudo-bulk RNA-seq samples (ED Fig. 1b)
```
Rscript Fig1b_tSNE_bulk.R
```

Check somatic mutations in this cohort (ED Fig. 1a)
```
Rscript SuppFig1a_oncoprint.R
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
Fig3e_CEBPA_ATAC.R
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
## Bulk RNA-seq quantification of B-ALL developmental states

We provide a series of notebooks outlining our approach to quantification of B-ALL developmental states in bulk RNA-seq data. This approach starts with: 
 - [Identification of differentially expressed genes across B-ALL developmental states](https://htmlpreview.github.io/?https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig5.0a_setup_BALL_Lineage_DE.nb.html)
 - [Identification of genes correlated with developmental state abundance across 85 pseudobulk profiles](https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig5.0b_NMFregression.R), filtered further by applying [adaptive thresholding within significantly correlated genes](https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig5.0d_MarkerGeneSelection.nb.html)
 - [Identification of differentially expressed genes across normal B cell development](https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig5.0c_BDev_pseudobulkDE.nb.html)

From the intersection of these genesets (differentially expressed across normal B cell development, differentially expressed across B-ALL cell developmental states, and correlated with developmental state abundance across pseudobulk profiles), we developed a regression based approach to [quantify B-ALL developmental states in bulk RNA-seq data](https://htmlpreview.github.io/?https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig5.0e_DevState_LASSOregression.nb.html). We found that for large scRNA-seq cohorts like this one (89 samples), this outperformed deconvolution with CIBERSORTx, BayesPrism, and DWLS for accurate quantification of various malignant cell states.

Notebooks and scripts related to methodology in Figure 5
```
Fig5.0a_setup_BALL_Lineage_DE.Rmd
Fig5.0b_NMFregression.R
Fig5.0c_BDev_pseudobulkDE.Rmd
Fig5.0d_MarkerGeneSelection.Rmd
Fig5.0e_DevState_LASSOregression.Rmd
```

==========================================================
## B-ALL Cohort Analysis and the B-ALL Multipotency Score 

Following successful quantification of B-ALL developmental states in bulk RNA-seq data, we applied this to a cohort of 2046 samples spanning pediatric and adult B-ALL. Principal component analysis revealed an axis of Multipotency vs Commitment, and we captured this axis through a 99-gene score which we termed the B-ALL Multipotency Score. The derivation and clinical characterization of this score are shown in this [R notebook](https://htmlpreview.github.io/?https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig5and6_MultipotencyScore_ClinicalAssociations.nb.html), with further clinical characterization done in this [jupyter notebook](https://htmlpreview.github.io/?https://github.com/gaoqs313/BALL-Single-Cell-Landscape/blob/main/Fig5and6_MultipotencyScore_ClinicalAssociations_extended.ipynb). These two notebooks cover the analyses for Figures 5-7 of the manuscript.

```
Fig5-7_MultipotencyScore_ClinicalAssociations.Rmd
Fig5-7_MultipotencyScore_ClinicalAssociations_extended.ipynb
```

