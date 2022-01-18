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

# STEP 1: Choosing and checking the currated miRNA-target database 
# Due to limatation of data avaliblty the miTarBase database is selcted for use.
# Human miTarBase interactions:
# hsa_MTI.xlsx saved in csv format 
# https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/download.php

# Checking the: file
setwd("C:/Users/modos/OneDrive - Norwich BioScience Institutes/Book Chapter - miRNA in IBD/scripts")
mitarbase_data <- read.csv("data/hsa_MTI.csv")
head(mitarbase_data)
dim(mitarbase_data)

# Based on this we 502 652 interactions. The IDs are ENTREZ Gene IDs and miRNA names which makes it easy to use.
length(unique(mitarbase_data$miRNA))
length(unique(mitarbase_data$Target.Gene))
# The database contains 2599 unique miRNA and 15 064 uniqe human genes.  




#STEP 2: Choosing and filtering the predicted miRNA-target database
# The following input fdiles are usedc for the pipline:
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
# website it is a trade off between sensitivity and specificty. We suggest using around a 
# 0.5 cut off for general miRNA predictions. The conserved site filtering and the PCT
# score should give a realtively good coverage. 

TS_filtered <- TargetScanData[TargetScanData$PCT > 0.5,]
dim(TS_filtered)

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
# The wheighted context+score tells the biochemical chance of certain miRNA-target gene site to be active. 
# If it is negative it is better. Here we can use directly the score or the percetiles function. We can chose 
# example the top 50% of miRNA -target gene peredctions or only certain genes. This is a trade off between 
# specificity and sesitivity like always. I have schosena  high specificity filter here. (Only the top 5% kept) 
TS_context_selected <- TargetScanSitePredictionsHuman[TargetScanSitePredictionsHuman$weighted.context...score.percentile>95, ]
plot3 <- qplot (TS_context_selected$weighted.context...score, 
                geom = "histogram",
                xlab = "Weighted context+ score",
                ylab = "Count",
                binwidth = 0.05)
plot3 + theme_classic() 
dim(TS_context_selected)
# This results 508 526 interactions. Note theese interactions are not miRNA - TG but miRNA 
# family target gene itneractions


# STEP 3: miRNA ID mapping 
TS_mirna_families <- read.csv("data/miR_Family_Info.txt", sep="\t")
head(TS_mirna_families)
#Collecting only huamn miRNAs. If you work in mouse or rats please change accordingly 
# checking the species ID here: https://www.ncbi.nlm.nih.gov/taxonomy 
TS_mirna_families_human <- TS_mirna_families[TS_mirna_families$Species.ID == 9606,]
#merging the two datasets. 
TS_merged <- merge(TS_filtered, TS_mirna_families_human, by.x = "miR.Family", by.y = "miR.family")
# This filtering results a large amount of interaction loss. It is bacause TS is based 
# on various species conservation infomration and not all mirNA family is present in humans.
dim(TS_merged)

TS_miRNANames = TS_merged$MiRBase.ID
version=checkMiRNAVersion(TS_miRNANames, verbose = TRUE)
# In the case of TargetScan the best version would be 21 

mitarbase_miRNANames <- mitarbase_data$miRNA
version=checkMiRNAVersion(mitarbase_miRNANames, verbose = TRUE)
# Luckily it is the same in the case of the miTarBase

#Now we can convert the two set of IDs to miRbase accesion numbers.
TS_miRNA_acession = miRNA_NameToAccession(TS_miRNANames,version = "v21")
TS_merged$Accession <- TS_miRNA_acession$Accession
mitarbase_miRNA_acession = miRNA_NameToAccession(mitarbase_miRNANames,version = "v21")
mitarbase_data$Accession <-  mitarbase_miRNA_acession$Accession
# We can add the representative tables the miRNA acessions. 

# STEP 4: Human ID mapping
#For simplicity we can map to ENTREZ GENE ID both dataset.
# ENSMEBLE mapping is hard due to version inconsistencies
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
entrez_gene_mapping <- AnnotationDbi::select(org.Hs.eg.db, keys=TS_merged$Gene.Symbol, columns='ENTREZID', keytype='SYMBOL')
TS_all <- merge(TS_merged, entrez_gene_mapping, by.x = "Gene.Symbol", by.y = "SYMBOL")

#STEP 5:  