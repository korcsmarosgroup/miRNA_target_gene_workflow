# Workflow code 2
# Analysing miRNA AND mRNA data
# Paper: https://pubmed.ncbi.nlm.nih.gov/30657881/  Gene and Mirna Regulatory Networks During Different Stages of Crohn's Disease
library("miRBaseConverter")
library("AnnotationDbi")
library("GEOquery")
library("affy")
library("miRLAB")
library("clusterProfiler")
library("igraph")
library("ggplot2")
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)
library("data.table")
library(limma)
library(umap)
library(tidyverse)
library(ReactomePA)
library(clusterProfiler)
library(ggpubr)

setwd("C:/Users/modos/OneDrive - Norwich BioScience Institutes/Book Chapter - miRNA in IBD/scripts/")

# Steps 1-4 Fromt he previous workflow.
# read combined TargetScan - mirtarbase prediction list from the previous workflow
mirna_TG <- read_csv("results/TS_mirtarbase.csv")

# Step 5
# Downloading expression data for the analyses
# First for mRNA
# load series and platform data from GEO
Sys.setenv("VROOM_CONNECTION_SIZE" = 199 * (1024^2)) # Setting max download to 200 MB
gset <- getGEO("GSE102133", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples X means groups excluded. 
# In this particular case we are comaparing control data with
# early crohn disease patinets.
# make sure they match between mRNA and miRNA!
gsms <- paste0(
  "000000000000XXX11111111111111XXXXXXXXXXXXXXXXXXXXX",
  "XXXXXXXXXXXXXXXXXXXXXXXXXXX"
)
sml <- strsplit(gsms, split = "")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[, sel]

# get groups
phenotypes <- gset@phenoData@data
control <- phenotypes %>%
  dplyr::filter(grepl("control", title)) %>%
  pull(geo_accession) # select first 4 so control and disease matches
control <- control[1:4]
lateCD <- phenotypes %>%
  dplyr::filter(grepl("late", title)) %>%
  pull(geo_accession)


# log2 transformation This is from GEO2R  Ithe script tries do decide 
# whether the downloaded data are log2 transformed.
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}

id2gn <- gset@featureData@data %>%
  dplyr::select(ID, Gene.symbol) %>%
  dplyr::filter(Gene.symbol != "")
id2gn$ID <- as.character(id2gn$ID)

# We prepeare table s for furtther worka nd correation calcualtions 
# between mRNA and miRNA levels.
mRNAdata <- ex %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(id2gn, by = c("ID" = "ID"))

mRNAdataControl <- mRNAdata %>%
  dplyr::select(control, Gene.symbol) %>%
  drop_na() %>%
  t()

mRNAdataCD <- mRNAdata %>%
  dplyr::select(lateCD, Gene.symbol) %>%
  drop_na() %>%
  t()

#Step 5 for miRNA 
# We have to redo it for miRNA
# It is the same study but the sample are used in all mammlian miRNA chip. 
# This will result a high amount of non-human miRNAs due to the hybridisation
# of sampe DNA to the probesets DNA on the CHIP.

gset <- getGEO("GSE102127", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL14613", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
# Note the samples are int he same order as in the mRNA and they are comming 
# form the smae patients. 
gsms <- paste0(
  "000000000000XXX11111111111111XXXXXXXXXXXXXXXXXXXXX",
  "XXXXXXXXXXXXXXXXXXXXXXXXXXX"
) # this needs to beckhecked that we have the same data.
sml <- strsplit(gsms, split = "")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[, sel]

# log2 transformation This we can leave out. We are using directly the genes et matrix.
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}

# get groups
phenotypes <- gset@phenoData@data
control <- phenotypes %>%
  dplyr::filter(grepl("control", title)) %>%
  pull(geo_accession) # select first 4
control <- control[1:4]
lateCD <- phenotypes %>%
  dplyr::filter(grepl("late", title)) %>%
  pull(geo_accession)


# Filter data
exFiltered <- ex %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  dplyr::filter(grepl("hsa-", ID))

miRNAdataControl <- exFiltered %>%
  dplyr::select(control, ID) %>%
  drop_na() %>%
  t()

miRNAdataCD <- exFiltered %>%
  dplyr::select(lateCD, ID) %>%
  drop_na() %>%
  t()

# unify rownames
rownames(miRNAdataCD) <- c()
rownames(miRNAdataControl) <- c()
rownames(mRNAdataCD) <- c()
rownames(mRNAdataControl) <- c()
#Checking the size of the data
dim(miRNAdataCD) # 2226
dim(mRNAdataCD) # 22195
dim(miRNAdataControl) # 2226
dim(mRNAdataControl) # 22195

# Step 8 Correaltion calculations with miRLAB

# join control-control and disease-disease tables
combinedControl <- cbind(mRNAdataControl, miRNAdataControl)
combinedCD <- cbind(mRNAdataCD[11:15, ], miRNAdataCD)

# move last row (ID) to first to be used as header by miRLAB
combinedControl <- combinedControl[c(5, 1:4), ]
combinedCD <- combinedCD[c(5, 1:4), ]



# write to file to be read by miRLAB
# the package will not read dataframes, only from file
write.table(combinedControl, "results/control_mrna_mirna.csv", sep = ",", col.names = F, row.names = F, quote = F)
write.table(combinedCD, "results/CD_mrna_mirna.csv", sep = ",", col.names = F, row.names = F, quote = F)

# calcaulting miRNA- target gene connections with miRLAB Pearson correation function
# The rows telling which are miRNA and which are mRNA. 
controlResults <- Pearson(datacsv = "results/control_mrna_mirna.csv", 1:22195, 22196:24421)
CDresults <- Pearson(datacsv = "results/CD_mrna_mirna.csv", 1:22195, 22196:24421)

# let's look at the distribution of correlation results
hist(controlResults)
hist(CDresults)

# keeping the bottom quartile of anticorrelation data, between -1 and -0.75
# converting it to dataframe first
# unifying IDs by changing . to - and similar

controlResultsFormatted <- controlResults %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(key = "key", value = "value", -rowname) %>%
  dplyr::filter(value <= -0.75) %>%
  mutate(rowname = str_replace_all(rowname, "[.]", "-")) %>%
  mutate(rowname = str_remove(rowname, "star")) %>%
  mutate(rowname = str_remove(rowname, "_st")) %>%
  mutate(rowname = str_remove(rowname, "_x"))


CDresultsFormatted <- CDresults %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  gather(key = "key", value = "value", -rowname) %>%
  dplyr::filter(value <= -0.75) %>%
  mutate(rowname = str_replace_all(rowname, "[.]", "-")) %>%
  mutate(rowname = str_remove(rowname, "star")) %>%
  mutate(rowname = str_remove(rowname, "_st")) %>%
  mutate(rowname = str_remove(rowname, "_x"))


# filter results using the TargetScan - mirtarbase data
controlResultsFiltered <- left_join(mirna_TG, controlResultsFormatted, by = c("miRNA" = "rowname", "Gene.Symbol" = "key")) %>% drop_na()
CDResultsFiltered <- left_join(mirna_TG, CDresultsFormatted, by = c("miRNA" = "rowname", "Gene.Symbol" = "key")) %>% drop_na()
write_csv(controlResultsFiltered, "results/control_edgelist.csv")
write_csv(CDResultsFiltered, "results/CD_edgelist.csv")


# Step 9 - Network visualisation and topological analysis
# Looking at the intersection and differences of the control and CD interactions

matchingInteractions <- intersect(controlResultsFiltered[, 1:2], CDResultsFiltered[, 1:2]) # 1243 common interactions
controlSpecific <- setdiff(controlResultsFiltered[, 1:2], CDResultsFiltered[, 1:2]) # 8744 interactions specific to control
cdSpecific <- setdiff(CDResultsFiltered[, 1:2], controlResultsFiltered[, 1:2]) # 8254 interactions specific to CD

# Get top 10 highest degree nodes from the specific networks, and run enrichment analysis on them

top10nodesCD <- cdSpecific %>%
  group_by(Gene.Symbol) %>%
  tally() %>%
  arrange(desc(n)) %>%
  head(10) %>%
  pull(Gene.Symbol)


top10nodesControl <- controlSpecific %>%
  group_by(Gene.Symbol) %>%
  tally() %>%
  arrange(desc(n)) %>%
  head(10) %>%
  pull(Gene.Symbol)

#Step 11 Gene ontology enrichmnet for the specific networks and hubs

# translate gene symbols to entrez
top10nodesCDentrez <- bitr(top10nodesCD, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
top10nodesControlEntrez <- bitr(top10nodesControl, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
background <- bitr(unique(mirna_TG$Gene.Symbol), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# Reactome
cdReactome <- enrichPathway(gene = top10nodesCDentrez$ENTREZID, pvalueCutoff = 0.05, readable = T, universe = background$ENTREZID)
cdReactomeResults <- cdReactome@result
controlReactome <- enrichPathway(gene = top10nodesControlEntrez$ENTREZID, pvalueCutoff = 0.05, readable = T, universe = background$ENTREZID)
controlReactomeResults <- controlReactome@result

# KEGG
cdKEGG <- enrichKEGG(gene = top10nodesCDentrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, universe = background$ENTREZID)
controlKEGG <- enrichKEGG(gene = top10nodesControlEntrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, universe = background$ENTREZID)


# GO
cdGO <- enrichGO(gene = top10nodesCDentrez$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.05, universe = background$ENTREZID)
controlGO <- enrichGO(gene = top10nodesControlEntrez$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.05, universe = background$ENTREZID)

# plots for CD
cd1 <- dotplot(cdReactome) # no enrichment in Reactome
cd2 <- dotplot(cdKEGG)
cd3 <- dotplot(cdGO)

fig1 <- ggarrange(cd2, cd3, ncol = 2, labels = c("KEGG", "GO"))

# plots for control
con1 <- dotplot(controlReactome) # no enrichment in Reactome
con2 <- dotplot(controlKEGG)
con3 <- dotplot(controlGO)

fig2 <- ggarrange(con3,
  ggarrange(con1, con2, ncol = 2, labels = c("Reactome", "KEGG")),
  nrow = 2,
  labels = "GO"
)

# KEGG vs KEGG
fig3 <- ggarrange(cd2, con2, ncol = 2, labels = c("CD", "Control"))

# GO vs GO
fig4 <- ggarrange(cd3, con3, ncol = 2, labels = c("CD", "Control"))

#A simmilar anlaysis can be used for the specific networks. 
# We put in the CD specific network here.
cd_specifc_network_targets <- cdSpecific$Gene.Symbol
cd_specifc_network_targets_translated <- bitr(cd_specifc_network_targets, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

#Reactome
cdReactomeall <- enrichPathway(gene = cd_specifc_network_targets_translated$ENTREZID, pvalueCutoff = 0.05, readable = T, universe = background$ENTREZID)
cdReactomeallResults <- cdReactomeall@result
cd1 <- dotplot(cdReactomeall) 
cd1

#KEGG
cdKEGGall <- enrichKEGG(gene = cd_specifc_network_targets_translated$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, universe = background$ENTREZID)
cd2 <- dotplot(cdKEGGall)
cd2
#Gene Ontology
cdGOall <- enrichGO(gene = cd_specifc_network_targets_translated$ENTREZID, 
                    OrgDb = "org.Hs.eg.db", ont = "BP", 
                    pvalueCutoff = 0.05, universe = background$ENTREZID)
cd3 <- dotplot(cdGOall)
cd3

fig5 <-  ggarrange(cd1, cd2, cd3, 
                   ncol = 2,nrow = 2, 
                   labels = c("Reactome", "KEGG", "Gene Ontology Biological Process"))
fig5
