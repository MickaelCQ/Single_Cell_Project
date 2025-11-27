library(dplyr)
library(Seurat)
# parallel computation and memory allocation for Seurat
library(future)
library(ComplexHeatmap)
plan(sequential)
options(future.globals.maxSize=10*1024**3)
p5 <- read.csv("/home/DROPRET/Téléchargements/CAFs/CAFs/GSM4805570_CountsMatrix_20G00953M_TN.txt.gz",
               sep="\t")
p4a <- read.csv("/home/DROPRET/Téléchargements/CAFs/CAFs/GSM4805566_CountsMatrix_19G02977A_TN.txt.gz",
                sep="\t")
p4b <- read.csv("/home/DROPRET/Téléchargements/CAFs/CAFs/GSM4805568_CountsMatrix_19G02977B_TN.txt.gz",
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
head(caf)





library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)



# high mitochondrial gene % is likely to indicate dead cells
# extract from ensdb
ensdb.genes <- genes(EnsDb.Hsapiens.v86)
# mito gene from MT chromosome
MT.names <- ensdb.genes[seqnames(ensdb.genes) == "MT"]$gene_name

# get the count matrix
counts <- GetAssayData(caf, "counts")

# data for cleaning data next
umi.tot   <- colSums(counts)
gene.tot  <- colSums(counts > 0)

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

# cleaning
counts <- counts[, !bad]


good.genes <- rowSums(counts > 1) >= 0.01 * ncol(counts)
counts <- counts[good.genes, ]
dim(counts)

caf <- CreateSeuratObject(counts = counts,
                          project = "CAF_clean",
                          min.cells = 0,
                          min.features = 0)
caf

