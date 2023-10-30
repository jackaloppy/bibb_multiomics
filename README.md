# bibb_multiomics
This is a repo for analyzing Dr.Bibb's multiomics project


### Data structures
- 4 scRNA sequencing data (C1,C2,T1,T2)
- 4 scATAC sequencing data (C1,C2,T1,T2)

FASTQ files are stored on Box foder. According to the SOP, they are sequencing data produced from 10X Genomics Chromium Single Cell Multiome ATAC + Gene Expression. 


### Preprocessing steps
I first ran [Cell Ranger ARC pipeline](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc) (v2.0.2) against mouse reference (mm10 Reference - 2020-A-2.0.0) provided by 10x. I paired the scRNA and scATAC up for each sample (i.e. C1-scATAC with C1-scRNA) and used the _cellranger_arc_template.slurm_ in [folder](/scripts/preprocessing/). The QC reports (_c1_arc.html_ in [folder](https://jackaloppy.github.io/bibb_multiomics/qc_reports/c1_arc_multiome.html)) showed very poor results. 

This kept me thinking that maybe the C1-scATAC and C1-scRNA (and other pairs) were not generated from the same multiome experiment. 10X has the option to analyze only the ATAC data  or only the Gene Expression data from the single cell multiome experiment (See [1](https://kb.10xgenomics.com/hc/en-us/articles/360061165691-Can-I-analyze-only-the-ATAC-data-from-my-single-cell-multiome-experiment-) and [2](https://kb.10xgenomics.com/hc/en-us/articles/360059656912-Can-I-analyze-only-the-Gene-Expression-data-from-my-single-cell-multiome-experiment-)). So I analyzed them separately according to the instructions from 10X in HPC (see _c1_scatac.slurm_ and _c1_scrna.slurm_ in [folder](/scripts/preprocessing/)). However, both of them yielded poor results based on the QC reports (see _c1_scatac_multiome.html_ and _c1_scrna_multiome.html_ in [folder](/qc_reports)). 

The first thing come up to me is the reference species. Using mouse reference is mentioned in the SOP on the CyVerse folder, even though they are rat samples. 10x only provides pre-built humand and mouse references but it has the [option](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr) to build a custom reference, or we can find one online (though it's hard to tell the quality). 

The QC reports contain UMAP automatically generated from the pipeline, though they looks pretty noisy. I will analyze the pipeline outputs in Seurat to see if I can further clean the data and produce some good UMAPs. But I think we should address the low mapping rate first.
