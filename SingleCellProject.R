library(dplyr)
library(Seurat)
# parallel computation and memory allocation for Seurat
library(future)
library(ComplexHeatmap)
plan(sequential)
options(future.globals.maxSize=10*1024**3)

p5 <-read.csv("/home/DROPRET/Documents/git/Single_Cell_Project/Data/GSM4805570_CountsMatrix_20G00953M_TN.txt.gz",
               sep="\t")
p4a <- read.csv("/home/DROPRET/Documents/git/Single_Cell_Project/Data/GSM4805566_CountsMatrix_19G02977A_TN.txt.gz",
                sep="\t")
p4b <- read.csv("/home/DROPRET/Documents/git/Single_Cell_Project/Data/GSM4805568_CountsMatrix_19G02977B_TN.txt.gz",
                sep="\t")

caf.data <- data.matrix(cbind(p5,p4a,p4b))
ens <- read.csv("~/Téléchargements/CAFs/ensembl-38-108-genes.txt", sep="\t")
ens2symb <- setNames(ens$Gene.name, ens$Gene.stable.ID)
ens2type <- setNames(ens$Gene.type, ens$Gene.stable.ID)
symbols <- ens2symb[rownames(caf.data)]
types <- ens2type[rownames(caf.data)]
good <- types=="protein_coding" & !is.na(symbols) & !duplicated(symbols)
sum(good)
symb <- symbols[good]
caf.data <- caf.data[good,]
rownames(caf.data) <- symb
caf <- CreateSeuratObject(counts=caf.data, project="cafs",
                          min.cells=0.01*ncol(caf.data), min.features=1000)
caf
# 12561 features across 3773 samples within 1 assay 
# Active assay: RNA (12561 features, 0 variable features)





library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)


# high mitochondrial gene % is likely to indicate dead cells
# extract from ensdb
ensdb.genes <- genes(EnsDb.Hsapiens.v86)
# mito gene from MT chromosome
MT.names <- ensdb.genes[seqnames(ensdb.genes) == "MT"]$gene_name

# get the count matrix
counts <- GetAssayData(caf, "RNA")

# data for cleaning data next
umi.tot   <- colSums(counts)
gene.tot  <- colSums(counts > 0)

sum(rownames(counts) %in% MT.names)

# percent of mito
mito.pc <- colSums(counts[rownames(counts) %in% MT.names, ]) / umi.tot * 100
hist(mito.pc, breaks=50)

# clean if lore than 50000 UMIs
bad.high <- umi.tot > 50000
# clean if less than 1000 distinct gene
bad.low  <- gene.tot < 1000
# clean if more than 50% of mito
bad.mito <- mito.pc > 50

bad <- bad.high | bad.low | bad.mito

' # to check the quantity of removed element
table(bad.high)
table(bad.low)
table(bad.mito)
table(bad)
'
# cleaning
# Keep cells with at least 1000 distinct genes detected in 1% of the cells
counts <- counts[, !bad.low]
good.genes <- rowSums(counts > 1) >= 0.01 * ncol(counts)
counts <- counts[good.genes, ]

# Then eliminate cells with >50,000 total UMIs and >50% mitochondrial genes.
bad <- bad[!bad.low]         # resize bad to match remaining cells
counts <- counts[, !bad]

dim(counts) #[1] 8958 3728

caf <- CreateSeuratObject(counts = counts,
                          project = "CAF_clean",
                          min.cells = 0,
                          min.features = 0)
caf
# 8958 features across 3728 samples within 1 assay 
# Active assay: RNA (8958 features, 0 variable features)
# => features = gene, samples = cell. Changing the data before remove cell, but add gene (or reverso). its from line 
# good.genes <- rowSums(counts > 1) >= 0.01 * ncol(counts).


# elementary normalization
# two choice : passing by a matrix or by a seurat object. Seurat object can be really faster for pca than a matrix (5 sec vs 10 min). 

# Calcul du pourcentage mitochondrial pour régression
caf[["percent.mt"]] <- PercentageFeatureSet(
  caf, pattern = "^MT-"
)

# Normalisation moderne SCTransform
ncounts <- SCTransform(
  caf, # applique SCTransform
  vars.to.regress = "percent.mt", # régresse l'effet des gènes mitochondriaux
  verbose = TRUE # affiche la progression
)



'norm <- log2(1 + sweep(counts, 2, colSums(counts), "/")*1e5)
ncounts <- CreateSeuratObject(counts = counts)
ncounts[["norm"]] <- CreateAssayObject(counts = norm)
DefaultAssay(ncounts) <- "norm"
ncounts <- FindVariableFeatures(ncounts, selection.method = "vst", nfeatures = 2000)
ncounts <- ScaleData(ncounts)
'


# realizing the pca for the data
ncounts <- RunPCA(ncounts, npcs = 100, features = VariableFeatures(ncounts))
DimPlot(ncounts, reduction = "pca")
ElbowPlot(ncounts, ndims = 100)

# pca is needed for realising the tsne
ncountsTSNE <- RunTSNE(ncounts, dims=1:10 , perplexity=30)
DimPlot(ncountsTSNE, reduction = "tsne" , label = FALSE )


ncountsUMAP <- RunUMAP(ncounts, dims=1:10)
DimPlot(ncountsUMAP, reduction = "umap" , label=TRUE )


ncountsUMAP<-FindNeighbors(ncountsUMAP,dims=1:10)
ncountsUMAP <- FindClusters ( ncountsUMAP , resolution =0.15)



ncountsUMAP.markers <- FindAllMarkers ( ncountsUMAP , only.pos = TRUE , min.pct =0.25 ,
                                   logfc.threshold =0.25)
ncountsUMAP.markers %>%
  group_by ( cluster ) %>%
  slice_max ( n =2 , order_by = avg_log2FC )
DimPlot(ncountsUMAP, reduction = "umap" , label=TRUE )


signatureGene <- c(
  "PLN","SORBS2","PHLDA2","SNCG","MT1M","MYH11",
  "PTGDS","FBLN1","DCN","LUM","COL1A1","LTBP2",
  "FABP5","HIGD1B","AGT","RGS5","CPE","SSTR2"
)

signatureGeneMarker <- ncountsUMAP.markers %>%
  dplyr::filter(gene %in% signatureGene)
signatureGeneMarker$gene <- factor(signatureGeneMarker$gene, levels = signatureGene)

signatureGeneMarker <- signatureGeneMarker %>%
  arrange(gene)
DoHeatmap(ncountsUMAP, features = signatureGene) + NoLegend()


DimPlot(ncountsUMAP, reduction = "umap" , label=TRUE )
FeaturePlot (ncountsUMAP, features = c ("PLN","SORBS2","PHLDA2","SNCG","MT1M","MYH11")) # 6 for each gene
FeaturePlot (ncountsUMAP, features = c ("PTGDS","FBLN1","DCN","LUM","COL1A1","LTBP2")) # 6 for each gene
FeaturePlot (ncountsUMAP, features = c ("FABP5","HIGD1B","AGT","RGS5","CPE","SSTR2")) # 6 for each gene























