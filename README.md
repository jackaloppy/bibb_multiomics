# bibb_multiomics
This is a repo for analyzing Dr.Bibb's multiomics project


### Data structures
- 4 scRNA sequencing data (C1,C2,T1,T2)
- 4 scATAC sequencing data (C1,C2,T1,T2)
- Not sure if those are paired sample or not, but the QC reports for scATAC were not ideal, so I analyzed them separately.

FASTQ files are stored on Box foder. According to the SOP, they are sequencing data produced from 10X Genomics Chromium Single Cell Multiome ATAC + Gene Expression on rat. 

### Rat reference libraries construction and preprocessing
Raw data were initially ran against mouse reference as suggested by the SOP, but the mapping rates were extremely poor for all the samples. We first checked the FastQC reports to assess the quality of FASTQ files:
  - FastQC reports
    - C1-scRNA: [R1](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1_S52_R1_001_fastqc.html) and [R2](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1_S52_R2_001_fastqc.html)
    - C1-scATAC: [I1](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1-ATAC_S1_L001_I1_001_fastqc.html), [R1](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1-ATAC_S1_L001_R1_001_fastqc.html), [R2](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1-ATAC_S1_L001_R2_001_fastqc.html), and [R3](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1-ATAC_S1_L001_R3_001_fastqc.html)
From the FastQC reports we can see that the quality of sequencing is good, so we need to seek good reference transcriptome.

Following the [tutorial](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr) from 10x, I built the rat reference libraries based on *mRatBN7.2* assembly. I'm not sure if those are paired samples, so I first analyzed the scRNA-seq and scATAC-seq separately. The QC reporst are below:

- Cell Ranger scRNA-seq
  - C1: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c1_scrna.html)
  - C2: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c2_scrna.html)
  - T1: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t1_scrna.html)
  - T2: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t2_scrna.html)
- Cell Ranger scATAC-seq
  - C1: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c1_scatac.html)
  - C2: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c2_scatac.html)
  - T1: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t1_scatac.html)
  - T2: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t2_scatac.html)

From the QC reports, we can see that the quality of scRNA-seq is pretty good, but the scATAC-seq have some warnings. So I think it's better to analyze separately, so that the scATAC-seq won't affect the good quality of the scRNA-seq results. Nonetheless, I also did the Cell Ranger ARC (GEX + ATAC) pipeline, the QC reports are shown below:

- Cell Ranger ARC (GEX + ATAC):
  - C1: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c1_arc.html)
  - C2: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c2_arc.html)
  - T1: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t1_arc.html)
  - T2: [html](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t2_arc.html)


### scRNA analysis
Next I loaded the scRNA-seq outputs into R and used Seurat to generate the UMAP. (See [script](/scripts/analysis/scrna_analysis.R)). After integrating the four samples, I did some visualizations on data quality:

![scRNA Number of cells per sample](/results/ncells_scrna.png)
![Number of genes per cell](/results/ngenepercell_scrna.png)
![scRNA Number of cells per sample](/results/ncells_scrna.png)

See more in [folder](/results). Then I ran UMAP and annotated the clusters based on the following gene markers:

- Neurons: Gad1, Gad2, Grin2a
- Astrocytes: Gja1, Aqp4
- Oligodendrocytes: Plp1, Mag
- Microglia: Arhgap15, Csf1r
- OPC: Vcan, Pdgfra

Here is the annotated UMAP:
![scRNA annotated UMAP](/results/annotated_umap_scrna.png)

### scATAC analysis
I loaded the scATAC-seq outputs into R and used Signac to generate the UMAP. (See [script](/scripts/analysis/scatac_analysis.R)). I had problems annotating genes to ranges for the data, so I'm only able to integrate and generate UMAP wihtout gene annotation information:

![scATAC UMAP](/results/umap_scatac.png)

The loupe files generated by 10X Cell Ranger ATAC pipeline are able to show UMAP with gene annotation information. I'm currently working on aggregating the ATAC outpus using Cell Ranger ATAC aggr pipeline. Then I can work on the loupe file on Loupe Broswer to get the annotated ATAC UMAP.


