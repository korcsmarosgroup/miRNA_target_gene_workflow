# STEP 0: PInstallation and loading the necesearry packages.
# Please call in the isntall R package to isntall all the necesearry packages.
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



# STEP 1: Choosing and checking the currated miRNA-target database 
# Due to limatation of data avaliblty the miTarBase database is selcted for use.
# Human miTarBase interactions:
# hsa_MTI.xlsx saved in csv format 
# https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/download.php

# Checking the: file
setwd("~/OneDrive - Norwich BioScience Institutes/Documents/mirna_genes/code/")
mitarbase_data <- read.csv("data/hsa_MTI.csv")
head(mitarbase_data)
dim(mitarbase_data)

# Based on this we 502 652 interactions. The IDs are ENTREZ Gene IDs and miRNA names which makes it easy to use.
length(unique(mitarbase_data$miRNA))
length(unique(mitarbase_data$Target.Gene))
# The database contains 2599 unique miRNA and 15 064 uniqe human genes.  

#STEP 2: Choosing and filtering the predicted miRNA-target database
# The following input files are used for the pipline:
# Prediction based miRNA - target gene data:
# TargetScan 5.2 human miRNA family to target gene interactions: Predicted_Targets_Info.default_predictions.txt 
# TargetScan family to mirBase ID converting table: miR_Family_Info.txt
# Donwloaded: http://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi
# miRWalk interactions: hsa_miRWalk_3UTR.txt, hsa_miRWalk_5UTR.txt, hsa_miRWalk_CDS.txt
# Downloaded from: http://mirwalk.umm.uni-heidelberg.de/resources/
TargetScanData <- read.csv("data/Predicted_Targets_Info.default_predictions.txt", sep="\t")
dim(TargetScanData)
head(TargetScanData)
# TargetScan miR family interactions contain 718 223 interactions.
# We can check the PCT scores. First converting them to numeric values.
TargetScanData$PCT <-as.numeric(TargetScanData$PCT)
plot <- qplot (TargetScanData$PCT, 
               geom = "histogram",
               xlab = "Target Scan PCT score",
               ylab = "Count",
               binwidth = 0.05)
plot + theme_classic() 
#After looking at the PCT scores users need to choose a cut-off. Based on the TargetScan 
# website it is a trade off between sensitivity and specificity. We suggest using around a 
# 0.5 cut off for general miRNA predictions. The conserved site filtering and the PCT
# score should give a realtively good coverage. 

TS_filtered_pct <- TargetScanData[TargetScanData$PCT > 0.5,]
dim(TS_filtered_pct)

# TargetScan has also the Biochemistry based predictions. Here the instead of the ctext plus scores are used. 
# This way every prediction is a predicted site.
TargetScanSitePredictions<- read.csv("data/Predicted_Targets_Context_Scores.default_predictions.txt", sep="\t")
dim(TargetScanSitePredictions)
head(TargetScanSitePredictions)
TargetScanSitePredictions$weighted.context...score <- as.numeric(TargetScanSitePredictions$weighted.context...score)
# We have many non-human predictions which we can exclude
TargetScanSitePredictionsHuman <- TargetScanSitePredictions[TargetScanSitePredictions$Gene.Tax.ID==9606,]
dim(TargetScanSitePredictionsHuman)
plot2 <- qplot (TargetScanSitePredictionsHuman$weighted.context...score, 
               geom = "histogram",
               xlab = "Weighted context+ score",
               ylab = "Count",
               binwidth = 0.05)
plot2 + theme_classic() 
# The wheighted context+ score tells the biochemical chance of certain miRNA-target gene site to be active. 
# If it is negative it is better. Here we can use directly the score or the percetiles function. We can chose 
# example the top 50% of miRNA -target gene peredctions or only certain genes. This is a trade off between 
# specificity and sensitivity like always. I have chosen a  high specificity filter here. (Only the top 5% kept) 
TS_filtered_context <- TargetScanSitePredictionsHuman[TargetScanSitePredictionsHuman$weighted.context...score.percentile>95, ]
plot3 <- qplot (TS_filtered_context$weighted.context...score, 
                geom = "histogram",
                xlab = "Weighted context+ score",
                ylab = "Count",
                binwidth = 0.05)
plot3 + theme_classic() 
dim(TS_filtered_context)
# This results 48 365 interactions. Note these interactions are miRNA - TG interactions,
# comapred to the pct interactions. 


# STEP 3: miRNA ID mapping 
#First we work with the PCT part of the TargetScan
TS_mirna_families <- read.csv("data/miR_Family_Info.txt", sep="\t")
head(TS_mirna_families)
#Collecting only huamn miRNAs. If you work in mouse or rats please change accordingly 
# checking the species ID here: https://www.ncbi.nlm.nih.gov/taxonomy 
TS_mirna_families_human <- TS_mirna_families[TS_mirna_families$Species.ID == 9606,]
#merging the two datasets. 
TS_merged_pct <- merge(TS_filtered_pct, TS_mirna_families_human, by.x = "miR.Family", by.y = "miR.family")
# This filtering results a large amount of interaction loss. It is because TS is based 
# on various species conservation infomration and not all mirNA family is present in humans.
#The next step we can concatanate the two data frame. 
miRNA_tg_pct <- data.frame("MiRBase.ID"  = TS_merged_pct$MiRBase.ID, "Gene.Symbol" = TS_merged_pct$Gene.Symbol)
dim(miRNA_tg_pct)
miRNA_tg_context <- data.frame("MiRBase.ID" = TS_filtered_context$miRNA, "Gene.Symbol" = TS_filtered_context$Gene.Symbol)
dim(miRNA_tg_context)

#The next step can have multiple solution. Depending of the researcher's aim 
# and their aim of specificity versus sensitivity. If we want to have a 
# specific predicted dataset we can use only those interactions wqhich are in 
# both dataframe. If our aim is to have high snesitivity then we can use the 
# union of the two datasets. 
TS_merged <- intersect(miRNA_tg_pct[, 1:2], miRNA_tg_context[, 1:2])



TS_miRNANames = TS_merged$MiRBase.ID
version=checkMiRNAVersion(TS_miRNANames, verbose = TRUE)
# In the case of TargetScan the best version would be 21 

mitarbase_miRNANames <- mitarbase_data$miRNA
version=checkMiRNAVersion(mitarbase_miRNANames, verbose = TRUE)
# Luckily it is the same in the case of the miTarBase

#Now we can convert the two set of IDs to miRbase accesion numbers.
TS_miRNA_acession = miRNA_NameToAccession(TS_miRNANames,version = "v21")
TS_merged$Accession <- TS_miRNA_acession$Accession
TS_miRNA_acession <- miRNA_NameToAccession(TS_context_selected$miRNA, version = "v21")
TS_context_selected$Accession <- TS_miRNA_acession$Accession
mitarbase_miRNA_acession = miRNA_NameToAccession(mitarbase_miRNANames,version = "v21")
mitarbase_data$Accession <-  mitarbase_miRNA_acession$Accession
# We can add the representative tables the miRNA acessions. 

# STEP 4: Human ID mapping
#For simplicity we can map to ENTREZ GENE ID both dataset.
# ENSMEBLE mapping is hard due to version inconsistencies
entrez_gene_mapping <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(TS_merged$Gene.Symbol), columns='ENTREZID', keytype='SYMBOL')
entrez_gene_mapping_TS <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(TS_context_selected$Gene.Symbol), columns='ENTREZID', keytype='SYMBOL')

# Converting the data.frame object to data.table objects for speed.
TS_merged <- data.table(TS_merged)
entrez_gene_mapping <- data.table(entrez_gene_mapping)
TS_all <- merge(TS_merged, entrez_gene_mapping, by.x = "Gene.Symbol", by.y = "SYMBOL", allow.cartesian=TRUE)
dim(TS_all)
head(TS_all, n=30)

# joining TS and mirtarbase data, writing to file, as it will be used in the next workflow
miRNA_TG1<-TS_merged %>% dplyr::select(MiRBase.ID, Gene.Symbol)
miRNA_TG2<-mitarbase_data %>% dplyr::select(miRNA, Target.Gene)
miRNA_TG <- rbind(miRNA_TG1, miRNA_TG2, use.names=F) %>% unique()
write_csv(miRNA_TG, 'results/TS_mirtarbase.csv')



# STEP 5 and STEP 6: Analysing the transcriptomic data and finding differentially expressed genes. 
# First the transcriptomic data without miRNA.
# load series and platform data from GEO
# The code is the standard GEO2R pipline. 
# The daset is comparing ileal mucosa of Crohn disease patinets versus control ielal mucosa.
# I have chosen only a part of the daset. Of course you can choose any other contrasts as well.
# Please see the paper here: https://pubmed.ncbi.nlm.nih.gov/30657881/
Sys.setenv("VROOM_CONNECTION_SIZE" = 199 * (1024^2)) #Setting max download to 200 MB
gset <- getGEO("GSE102133", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("000000000000XXX11111111111111XXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("control","late_stage_CD"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B",number=Inf)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="results/Differetially_expressed_control_vs_highUC.txt", row.names=F, sep="\t")

#The below are standard visulasiations of microarray transcriptmic experimnets. Feel free to exclude. 
# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

DEGs <- tT2[tT2$adj.P.Val<0.05 & (tT2$logFC>1 | tT2$logFC< -1),]

# Step 7: Common miRNAs targeting the DEGs.
# Here there are many things to use. e.g. only the manually currated mirNA -target
# database like here
mitarbase_Degs_mirNAs <- merge(DEGs, mitarbase_data, by.x= "Gene.ID", by.y= "Target.Gene..Entrez.Gene.ID.")
head(mitarbase_Degs_mirNAs)
dim(mitarbase_Degs_mirNAs)
#We can simplify the results
miRNATG_for_graph <- unique(mitarbase_Degs_mirNAs[c("miRNA", "Gene.symbol")])
dim(miRNATG_for_graph)

# We can use the joined TS data as well.
TS_Degs_mirNAs <- merge(DEGs, TS_all, by.x= "Gene.ID", by.y= "ENTREZID")
head(TS_Degs_mirNAs)
dim(TS_Degs_mirNAs)
#For celan viosualisation we use this network
miRNATG_TS_for_graph <- unique(TS_Degs_mirNAs[c("MiRBase.ID", "Gene.symbol")])
dim(miRNATG_TS_for_graph)

# Step 9: Visualising the network
# Network visulalisation with R is  hard. 
mirnas <-unique(miRNATG_TS_for_graph$MiRBase.ID)
genes <- unique(miRNATG_TS_for_graph$Gene.symbol)
length(mirnas)
length(genes)
nodenames <- c(mirnas,genes)
type_vector <- c(rep("miRNA", length(mirnas)),rep("gene", length(genes)))
type_df <- data.frame(Name = nodenames, Type = type_vector) #into this dF you can addtional attributes.

network <- graph_from_data_frame(miRNATG_TS_for_graph, vertices= type_df )
degree_data <- degree(
    network,
    v = V(network),
    mode = c("all"),
    loops = TRUE,
    normalized = FALSE
)
network <- set_vertex_attr(network, "degree", value=degree_data)
vertex_attr(network,"degree")
colrs <- c("green", "blue")
colrs <- setNames(colrs, c("miRNA", "gene"))
V(network)$color <- colrs[V(network)$Type]
V(network)$size <- V(network)$degree*0,5
plot(network,  
     edge.arrow.size=.1,
     edge.color="black", 
     vertex.frame.color="gray",
     edge.curved=0.2,
     vertex.label="",
     vertex.label.color="black",
     layout=layout_with_kk,
     width=0.01)

analytical_outcome <- data.frame(Name = vertex_attr(network, "name"), 
                                Type = vertex_attr(network, "Type"), 
                                Degree= vertex_attr(network, "degree"))

#At the end we can check which gfenesa re the hubs in the netework
top10_degs <- analytical_outcome %>%
    arrange(desc(Degree)) %>%    
    filter(Type == "gene") %>%
    head(10)
top10_mirnas <- analytical_outcome %>%
    arrange(desc(Degree)) %>%    
    filter(Type == "miRNA") %>%
    head(10)

write_graph(
    network,
    "results/TS_DEGS.ncol",
    format = "ncol")
