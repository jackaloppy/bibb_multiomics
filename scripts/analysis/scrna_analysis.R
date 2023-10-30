## SoupX (need to download whole /outs folder)
sc_m1 = load10X("./data/c1_scrna_filtered_feature_bc_matrix")


c1_scrna <- Read10X("./data/c1_scrna_filtered_feature_bc_matrix/")
c1_scrna <- CreateSeuratObject(counts = c1_scrna, project = "c1_scrna", min.cells = 3, min.features=200)

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
