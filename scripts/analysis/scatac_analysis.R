counts <- Read10X_h5("./data/t2_scatac/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "./data/t2_scatac/singlecell.csv",
  header = TRUE,
  row.names=1
)
assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":","-"),
  genome = "mRatBN7.2",
  fragments = "./data/t2_scatac/fragments.tsv.gz",
  min.cells = 1
)
t2 <- CreateSeuratObject(
  counts = assay,
  assay = "peaks",
  project = "T2",
  meta.data = metadata
)


## Add gene annotations to the object
library(GenomeInfoDb)
gtf <- rtracklayer::import("./data/Rattus_norvegicus.mRatBN7.2.110.filtered.gtf")
gene.coords <- gtf[gtf$type == 'gene']
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

Annotation(c1) <- gene.coords
Annotation(c2) <- gene.coords
Annotation(t1) <- gene.coords
Annotation(t2) <- gene.coords

# Preprocessing individually
t2 <- NucleosomeSignal(t2)
t2$nucleosome_group <- ifelse(t2$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(t2, region = '1-1-1000000')

t2 <- TSSEnrichment(t2, fast = FALSE)
t2$high.tss <- ifelse(t2$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(t2, group.by = 'high.tss') + NoLegend()

t2$pct_reads_in_peaks <- t2$peak_region_fragments / t2$passed_filters * 100
t2$blacklist_ratio <- t2$blacklist_region_fragments / t2$peak_region_fragments


VlnPlot(
  object = merged,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 
               'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5,
  split.by = "dataset"
)

c1 <- FindTopFeatures(c1, min.cutoff = 10)
c1 <- RunTFIDF(c1)
c1 <- RunSVD(c1)


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
integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3)
DimPlot(integrated, split.by="dataset", group.by = "seurat_clusters")

gene.activities <- GeneActivity(integrated)

integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)

DefaultAssay(integrated) <- 'RNA'
FeaturePlot(
  object = integrated,
  features = c('Gad2', 'Gja1', 'Plp1', 'Csf1r', 'Vcan'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

FeaturePlot(integrated, features = "Grin2a")
DimPlot(integrated, split.by="dataset")
HoverLocator(plot=DimPlot(integrated))
select.cells <- CellSelector(plot=DimPlot(integrated))
Idents(integrated, cells=select.cells) <- "Neurons"

integrated <- RenameIdents(integrated, '18' = "Neurons")
saveRDS(integrated, "./data/merged_scatac.rds")




