counts <- Read10X_h5("./data/c2_scatac/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "./data/c2_scatac/singlecell.csv",
  header = TRUE,
  row.names=1
)
c1_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":","-"),
  genome = "mRatBN7.2",
  fragments = "./data/c2_scatac/fragments.tsv.gz",
  min.cells = 1
)
c2 <- CreateSeuratObject(
  counts = c1_assay,
  assay = "peaks",
  project = "C2",
  meta.data = metadata
)




## Add gene annotations to the object
library(AnnotationHub)
library(dbplyr)
library(BiocFileCache)
ah = AnnotationHub()
rat <- query(ah, c("EnsDb","Rattus norvegicus", 110))
rat_ref <- rat[["AH113793"]]
library(GenomeInfoDb)
gtf <- rtracklayer::import("./data/Rattus_norvegicus.mRatBN7.2.110.filtered.gtf")
gene.coords <- gtf[gtf$type == 'gene']
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

annotations <- GetGRangesFromEnsDb(ensdb = rat_ref)
seqlevelsStyle(annotations) <- "NCBI"
genome(annotations) <- "mRatBN7.2"
Annotation(c1) <- annotations
Annotation(c2) <- annotations
Annotation(t1) <- annotations
Annotation(t2) <- annotations

# Preprocessing individually
c1 <- TSSEnrichment(c1, fast = FALSE)

c2$pct_reads_in_peaks <- c2$peak_region_fragments / c2$passed_filters * 100
c2$blacklist_ratio <- c2$blacklist_region_fragments / c2$peak_region_fragments
c1 <- TSSEnrichment(c1, fast=FALSE)

VlnPlot(
  object = c2,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

t2 <- FindTopFeatures(t2, min.cutoff = 10)
t2 <- RunTFIDF(t2)
t2 <- RunSVD(t2)


## Merge
c1$dataset <- "C1"
c2$dataset <- "C2"
t1$dataset <- "T1"
t2$dataset <- "T2"
merged <- merge(x = c1, y = list(c2, t1, t2), add.cell.ids = c("C1","C2","T1","T2")) 

## Preprocessing
merged <- FindTopFeatures(merged, min.cutoff = 10)
merged <- RunTFIDF(merged)
merged <- RunSVD(merged)
merged <- RunUMAP(merged, dims=2:30, reduction="lsi")
split <- SplitObject(merged, split.by="dataset")

## Integration
integration.anchors <- FindIntegrationAnchors(
  object.list = split,
  anchor.features = rownames(split),
  reduction = "rlsi",
  dims = 2:30
)

integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merged[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
DimPlot(integrated, group.by="dataset")
DimPlot(integrated, split.by="dataset")


Annotation(integrated) <- gene.coords
integrated <- TSSEnrichment(integrated, fast = FALSE)
gene.activities <- GeneActivity(integrated)
