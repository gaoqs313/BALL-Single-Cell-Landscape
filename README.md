# BALL-Single-Cell-Landscape

This is the repository for scripts used in our Genome Biology paper "Splicing-accessible coding 3’UTRs control protein stability and interaction networks"


==========================================================
## Prepare penultimate exon annotation

### Mouse Reference

Version: Mus_musculus.GRCm38.98 \
Script folder: 01_mouse_reference

Extract eligible transcript
```
perl extract_eligible_transcript.pl Mus_musculus.GRCm38.98.chr.gtf mm10.fa
cat StopCodonInLastExon.final.Mus_musculus.GRCm38.98.chr.gtf StopCodonInPenultimateExon.final.Mus_musculus.GRCm38.98.chr.gtf > final.Mus_musculus.GRCm38.98.chr.gtf
```

Calculate amino acid content
```
perl calculate_AA.pl
Rscript stat.R
```

Convert GTF to MISO GFF3 file
```
perl GenerateMisoSE.pl
```

Build MISO index
```
index_gff --index mm10.SE.gff3 indexed_SE_events/
```

### Human Reference

Version: Homo_sapiens.GRCh38.98 \
Script folder: 02_human_reference

Extract eligible transcript
```
perl extract_eligible_transcript.pl Homo_sapiens.GRCh38.98.chr.gtf hg38.fa
cat StopCodonInLastExon.final.Homo_sapiens.GRCh38.98.chr.gtf StopCodonInPenultimateExon.final.Homo_sapiens.GRCh38.98.chr.gtf > final.Homo_sapiens.GRCh38.98.chr.gtf
```

Calculate amino acid content
```
perl calculate_AA.pl
Rscript stat.R
```

Convert GTF to MISO GFF3 file
```
perl GenerateMisoSE.pl
```

Build MISO index
```
index_gff --index hg38.SE.gff3 indexed_SE_events/
```

### Liftover

Script folder: 03_liftover

```
bash run_liftOver.sh
```

==========================================================
## Process RNA-seq data

Script folder: 04_data_processing

### 1. download
```
bash run_download.sh
```

### 2. align using tophat2
```
bash tophat.sh
```

### 3. run miso
```
bash insert_len.sh
bash miso.sh
```
### 4. run kallisto
```
bash kallisto.sh
```

==========================================================
## Do statistics and plotting   

Script folder: 05_statistics

### Mouse_tissue

Figures: Fig 1B, 1C, 1E, S1B, S1C
```
perl combine_mouse_tissue_penultimate.pl
Rscript	mouse_tissue_penultimate_stat.R
```
 
Figures: Fig S1F
```
Rscript combine_plot.R
```

### Mouse_neuron

Figures: Fig 1D, S1D, S1E
```
Rscript neuron_stat.R
```

Figures: Fig S1G, S1H
```
bash run_motif.sh
Rscript gene_stat.R
```

### Human_tissue

Figures: Fig 4A, S4A, S4B
```
perl combine_human.pl
Rscript human_stat.R
```

### Note 
Fig 4B, S4D are schematic illustration using Inkscape \
Fig 4C, S4C are generated using `plot_overlap_v2.R` in 03_liftover \
Fig 3A, 4E are generated using `mouse_aa_stat.R` in 01_mouse_reference and `human_aa_stat.R` in 02_human_reference \
Fig 3D is a simple boxplot generated in 2015 ...
