library("miRBaseConverter")
library("AnnotationDbi")
library("GEOquery")
library("affy")
library("miRLAB")
library("clusterProfiler")
library("igraph")

# STEP 1:
# Due to limatation of data avaliblty the miTarBase database is selcted for use.
# Human miTarBase interactions:
# MTI.xlsx saved in csv format 
# https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/download.php

# Checking the:



# The following input fdiles are usedc for the pipline:
# Prediction based miRNA - target gene data:
# TargetScan 5.2 human miRNA family to target gene interactions: Predicted_Targets_Info.default_predictions.txt 
# TargetScan family to mirBase ID converting table: miR_Family_Info.txt
# Donwloaded: http://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi
# miRWalk interactions: hsa_miRWalk_3UTR.txt, hsa_miRWalk_5UTR.txt, hsa_miRWalk_CDS.txt
# Downloaded from: http://mirwalk.umm.uni-heidelberg.de/resources/


