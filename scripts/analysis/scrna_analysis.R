### SoupX (need to download whole /outs folder)
### Don't need to do it now. It's a refined process that we could do later to better clean up the dataset.
### Same as DoubletFinder
#test_c1 <- load10X("./data/c1_scrna_outs")
#test_c1 <- autoEstCont(test_c1)
#test_c1 <- adjustCounts(test_c1)
### No need to use SoupX because estimated global rho is 0.03, which is low.


## Integrative analysis in Seurat V5
options(Seurat.object.assay.version = "v5")

c1_scrna <- Read10X("./data/c1_scrna_filtered_feature_bc_matrix/")
c1_scrna <- CreateSeuratObject(c1_scrna,min.cells = 3, min.features = 200)
c1_scrna$sample <- "c1"
c2_scrna <- Read10X("./data/c2_scrna_filtered_feature_bc_matrix/")
c2_scrna <- CreateSeuratObject(c2_scrna,min.cells = 3, min.features = 200)
c2_scrna$sample <- "c2"
t1_scrna <- Read10X("./data/t1_scrna_filtered_feature_bc_matrix/")
t1_scrna <- CreateSeuratObject(t1_scrna,min.cells = 3, min.features = 200)
t1_scrna$sample <- "t1"
t2_scrna <- Read10X("./data/t2_scrna_filtered_feature_bc_matrix/")
t2_scrna <- CreateSeuratObject(t2_scrna,min.cells = 3, min.features = 200)
t2_scrna$sample <- "t2"

lupus.anchors <- FindIntegrationAnchors(object.list=list(c1_scrna,c2_scrna,t1_scrna,t2_scrna))
combined <- IntegrateData(anchorset=lupus.anchors)
DefaultAssay(combined) <- "integrated"


c1_scrna[["percent.mt"]] <- PercentageFeatureSet(c1_scrna, pattern = "^MT-")
VlnPlot(c1_scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(c1_scrna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


c1_scrna <- NormalizeData(c1_scrna)
c1_scrna <- FindVariableFeatures(c1_scrna, selection.method = "vst", nfeatures = 2000)
c1_scrna <- ScaleData(c1_scrna, features = rownames(c1_scrna))
c1_scrna <- RunPCA(c1_scrna, features = VariableFeatures(object = c1_scrna))

c1_scrna <- JackStraw(c1_scrna, num.replicate = 100)
c1_scrna <- ScoreJackStraw(c1_scrna, dims = 1:20)
ElbowPlot(c1_scrna)


c1_scrna <- FindNeighbors(c1_scrna, dims = 1:6)
c1_scrna <- FindClusters(c1_scrna, resolution = 0.2)
c1_scrna <- RunUMAP(c1_scrna, dims = 1:6)
DimPlot(c1_scrna, reduction = "umap")

FeaturePlot(c1_scrna, "Ctss")
