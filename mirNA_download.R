library("miRBaseConverter")
library("AnnotationDbi")
library("GEOquery")
library("affy")
library("miRLAB")
library("clusterProfiler")
library("igraph")
library("ggplot2")
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library("data.table")
library(limma)
library(umap)

gset <- getGEO("GSE102127", GSEMatrix =TRUE, AnnotGPL=FALSE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL14613", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("000000000000XXX11111111111111XXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXX") # this needs to beckhecked that we have the same data.
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation This we can leave out. We are using directly the genes et matrix.
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
