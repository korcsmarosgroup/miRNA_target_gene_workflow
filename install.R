# This is the installation notebook

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("miRBaseConverter")
BiocManager::install("AnnotationDbi")
BiocManager::install("GEOquery")
BiocManager::install("affy")
BiocManager::install("miRLAB")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("limma")
BiocManager::install("umap")

install.packages("igraph")
install.packages("ggplot2")
install.packages("data.table")