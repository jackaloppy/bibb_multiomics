### SoupX (need to download whole /outs folder)
### Don't need to do it now. It's a refined process that we could do later to better clean up the dataset.
### Same as DoubletFinder
# test_c1 <- load10X("./data/c1_outs/")
# test_c1 <- autoEstCont(test_c1)
# test_c1 <- adjustCounts(test_c1)
### No need to use SoupX because estimated global rho is 0.02, which is low.


## Integrative analysis in Seurat V5

c1_scrna <- Read10X("./data/c1_scrna_filtered_feature_bc_matrix/")
c1_scrna <- CreateSeuratObject(c1_scrna, project="C1", min.cells=3, min.feature=200)
#c1_scrna$sample <- "c1"
c2_scrna <- Read10X("./data/c2_scrna_filtered_feature_bc_matrix/")
c2_scrna <- CreateSeuratObject(c2_scrna, project="C2",min.cells=3, min.feature=200)
#c2_scrna$sample <- "c2"
t1_scrna <- Read10X("./data/t1_scrna_filtered_feature_bc_matrix/")
t1_scrna <- CreateSeuratObject(t1_scrna, project="T1",min.cells=3, min.feature=200)
#t1_scrna$sample <- "t1"
t2_scrna <- Read10X("./data/t2_scrna_filtered_feature_bc_matrix/")
t2_scrna <- CreateSeuratObject(t2_scrna, project="T2",min.cells=3, min.feature=200)
#t2_scrna$sample <- "t2"

merged_obj <- merge(x = c1_scrna, y=c(c2_scrna, t1_scrna, t2_scrna), add.cell.ids = c("C1","C2","T1","T2"))

merged_obj$log10GenesPerUMI <- log10(merged_obj$nFeature_RNA) / log10(merged_obj$nCount_RNA)
merged_obj$mitoRatio <- PercentageFeatureSet(merged_obj, pattern = "^Mt-")
merged_obj$mitoRatio <- merged_obj@meta.data$mitoRatio / 100
merged_obj$riboRatio <- PercentageFeatureSet(merged_obj, pattern="^Rp[sl]")
merged_obj$riboRatio <- merged_obj@meta.data$riboRatio / 100
merged_obj$hbRatio <- PercentageFeatureSet(merged_obj, pattern="Hb[^(p)]")
merged_obj$hbRatio <- merged_obj@meta.data$hbRatio / 100

metadata <- merged_obj@meta.data
metadata$cells <- rownames(metadata)

metadata$sample <- metadata$orig.ident
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells_scRNA")

options(scipen = 999)
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10(limits = c(100,NA)) + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+
  ggtitle("UMI counts (transcripts) per cell")

metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10(limits = c(200,NA)) + 
  ylab("Cell density") +
  geom_vline(xintercept = 300)+
  ggtitle("nGene per cell")

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  xlim(0.7, NA) +
  ggtitle("Gene detected per UMI (novelty score)")

metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) +
  ggtitle("Mitochondrial gene ratio")

metadata %>% 
  ggplot(aes(x=nUMI, y=nGene)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  facet_wrap(~sample) +
  ggtitle("UMIs vs. genes")

metadata %>% 
  ggplot(aes(color=sample, x=hbRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) +
  ggtitle("Mitochondrial gene ratio")


VlnPlot(merged_obj, features=c("mitoRatio","riboRatio","hbRatio"))

merged <- readRDS("./data/merged.rds")

## Integration
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
merged_obj <- IntegrateLayers(merged_obj, method=CCAIntegration, 
                              orig.reduction="pca", new.reduction="integrated.cca",
                              verbose=FALSE)

ElbowPlot(merged_obj)
merged_obj <- FindNeighbors(merged_obj, reduction="integrated.cca", dims = 1:14)
merged_obj <- FindClusters(merged_obj, resolution = 0.1, cluster.name="cca_clusters")
merged_obj <- RunUMAP(merged_obj, dims=1:14, reduction="integrated.cca", reduction.name="umap.cca")

DimPlot(merged_obj, reduction="umap.cca", group.by = "cca_clusters", split.by = "orig.ident")

FeaturePlot(merged_obj, reduction="umap.cca", features = "Cd79a")
DimPlot(merged_obj, reduction="umap.cca", split.by = "orig.ident")

HoverLocator(plot=DimPlot(merged_obj,reduction="umap.cca"))
select.cells <- CellSelector(plot=DimPlot(merged_obj,reduction="umap.cca"))
Idents(merged_obj, cells=select.cells) <- "OPC"

saveRDS(merged_obj, file="./data/merged.rds")
