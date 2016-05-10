#Preprocessing data first step

### RUV Analysis for Combined Chn+Australian (QLD) ALS samples ###
### Modified by Beben Benyamin from Maksimovic et al 2015 ####

#source("http://bioconductor.org/biocLite.R")

#biocLite("sva")
#biocLite("ROCR")
#biocLite("isva")
#biocLite("FlowSorted.Blood.450k")

library(missMethyl)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(sva)
library(ruv)
library(ROCR)
library(RColorBrewer)
library(matrixStats)
library(isva)
library(FlowSorted.Blood.450k)
source("~/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/utils.R")
source("~/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/functions.R")
ann.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
Xreact = read.csv(file="~/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/48639-non-specific-probes-Illumina450k.csv", 
                  stringsAsFactors=FALSE)

## Set directory
data="~/clusterdata/uqpwains/ibscratch/RAW-IDATS"
directory="~/clusterdata/uqpwains/ibscratch/QC-IDATS"

#Set to data directory
setwd(data)

### Load Combined target data
sample_file="SampleSheetPD1-10.csv.csv"

#twin data
#twin="~/ibscratch/MND/ExperimentalDesign/QLD/20141003_Wray/Analysis2/DuplicateTwin.txt"
#twin=read.table(twin, header=T)

## Specify covariates for which methylation data should be adjusted
covariates=c("predicted_sex","horvath_predicted_age","plate","array", 
             "array_pos","CD8T","CD4T","NK","Bcell","Mono","Gran")

## read in targets and idats
targets=read.csv(sample_file, header=T, stringsAsFactor=F)
targets$array <- as.factor(targets$array)
#targets=subset(targets, targets$Cohort == "QLD")

rgSet=load("RGset.Robject");
detP = detectionP(rgSet) ## calculate detection p-values

# Change directory to the analysis directory
setwd(directory) 
save(rgSet, file="RGset.RObject")

#######################################################################################

## Examine mean detection p-values for all samples
pdf("detection-p-values.pdf", width=14)
barplot(apply(detP,2,mean), main="mean detection p-values", col=as.factor(targets$Source), xaxt="none")
abline(h=0.01,col="red")
dev.off()

## Pre-process the data after excluding poor quality samples
set.seed(100)
mSetSw = preprocessSWAN(rgSet[,apply(detP,2,mean) < 0.01]) ## probe-type normalisation only
gmSetQ = preprocessQuantile(rgSet[,apply(detP,2,mean) < 0.01]) ## probe-type and btw array normalisation
gmSetQ = gmSetQ[match(featureNames(mSetSw),featureNames(gmSetQ)),] ## ensure probes are ordered the same

## Remove poor quality samples from targets info and detection p-values  
targets = targets[apply(detP,2,mean) < 0.01,]
detP = detP[,apply(detP,2,mean) < 0.01]

save(targets, file="targets.RObject")
save(detP, file="detP.RObject")
