# bibb_multiomics
This is a repo for analyzing Dr.Bibb's multiomics project


### Data structures
- 4 scRNA sequencing data (C1,C2,T1,T2)
- 4 scATAC sequencing data (C1,C2,T1,T2)

FASTQ files are stored on Box foder. According to the SOP, they are sequencing data produced from 10X Genomics Chromium Single Cell Multiome ATAC + Gene Expression. 

### QC reports
- Cell Ranger RNA (against mouse reference)
  - C1 : [link](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c1_scrna_web_summary.html)
  - C2 : [link](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c2_scrna_web_summary.html)
  - T1 : [link](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t1_scrna_web_summary.html)
  - T2 : [link](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t2_scrna_web_summary.html)
- Cell Ranger ATAC (against mouse reference)
  - C1 : [link](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c1_scatac_web_summary.html)
  - C2 : [link](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c2_scatac_web_summary.html)
  - T1 : [link](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t1_scatac_web_summary.html)
  - T2 : [link](https://jackaloppy.github.io/bibb_multiomics/qc_reports/t2_scatac_web_summary.html)
- Cell Ranger ARC (ATAC + RNA) (against mouse reference)
  - Only C1 was tested: [link](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c1_arc_multiome.html)

### Preprocessing steps
I first ran [Cell Ranger ARC pipeline](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc) (v2.0.2) against mouse reference (mm10 Reference - 2020-A-2.0.0) provided by 10x. I paired the scRNA and scATAC up for each sample (i.e. C1-scATAC with C1-scRNA) and used the _cellranger_arc_template.slurm_ in [folder](/scripts/preprocessing/). The QC reports (see _C1 Cell Ranger ARC (ATAC + RNA)_ above) showed very poor results. Sepcifically, the estimated cells are reaching the upper limit of 20k and doubles the estimated cells from RNA or ATAC pipeline, suggesting barcodes are not shared between RNA and ATAC pipelines. So they are probably individual (unpaired) samples.

10X has the option to analyze only the ATAC data or only the Gene Expression data from the single cell multiome experiment (See [1](https://kb.10xgenomics.com/hc/en-us/articles/360061165691-Can-I-analyze-only-the-ATAC-data-from-my-single-cell-multiome-experiment-) and [2](https://kb.10xgenomics.com/hc/en-us/articles/360059656912-Can-I-analyze-only-the-Gene-Expression-data-from-my-single-cell-multiome-experiment-)). So I analyzed them separately according to the instructions from 10X in HPC (see _c1_scatac.slurm_ and _c1_scrna.slurm_ in [folder](/scripts/preprocessing/)). However, all of them yielded poor mapping rates from the QC reports (see above). 

### next steps
- I will assess the FASTQ file qulaity score and map the histogram across all samples, because bad quality socre could result in poor hit, causing low mapping rate. I will also try to increase the mis-match threshold in the 10x Cell Ranger pipeline to see if that increases the mapping rate. If we showed that the low mapping rate is due to low quality score, then something bad happened during the sequencing steps and nothing we can remedy that.
  - FastQC reports
    - C1-scRNA: [R1](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1_S52_R1_001_fastqc.html) and [R2](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1_S52_R2_001_fastqc.html)
    - C1-scATAC: [I1](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1-atac_S1_L001_I1_001_fastqc.html), [R1](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1-atac_S1_L001_R1_001_fastqc.html), [R2](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1-atac_S1_L001_R2_001_fastqc.html), and [R3](https://jackaloppy.github.io/bibb_multiomics/FastQC_reports/C1-atac_S1_L001_R3_001_fastqc.html)
  - From the FastQC reports we can see that the quality of sequencing is good, so we need to seek good reference transcriptome.
- Another thing we can check is to use rat referecne. The SOP on the CyVerse folder suggests using the pre-built mouse reference from 10x, even though they are rat samples (Not sure the reasoning behind). However, 10x has the [option](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr) to build a custom reference. And there are rat reference library available oneline. I will see if I can incoporate that into the pipeline, to increase the mapping rate.

- I will also double check the pair matches. Specifically, check if C1(RNA) can matches C2(ATAC). 

- Rudra is checking the potential microbes contamination by alinging the fastq with microbiome.

### on hold
The QC reports contain UMAP automatically generated from the pipeline, though they looks pretty noisy. I will analyze the pipeline outputs in Seurat to see if I can further clean the data and produce some good UMAPs. But I think we should address the low mapping rate first.
