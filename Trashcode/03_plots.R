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
source("~/ibscratch/workdir/scripts/QCscript/utils.R")
source("~/ibscratch/workdir/scripts/QCscript/functions.R")
ann.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
Xreact = read.csv(file="~/ibscratch/workdir/scripts/QCscript/48639-non-specific-probes-Illumina450k.csv", 
                  stringsAsFactors=FALSE)

## Set directory
data="~/ibscratch/RAW-IDATS"
directory="~/ibscratch/QC-IDATS"

#Set to data directory
setwd(data)

### Load Combined target data
sample_file="SampleSheetPD1-10.csv"



#twin data
#twin="~/ibscratch/MND/ExperimentalDesign/QLD/20141003_Wray/Analysis2/DuplicateTwin.txt"
#twin=read.table(twin, header=T)

## Specify covariates for which methylation data should be adjusted
#covariates=c("predicted_sex","horvath_predicted_age","plate","array", 
#             "array_pos","CD8T","CD4T","NK","Bcell","Mono","Gran")

## read in targets and idats
targets=read.csv(sample_file, header=T, stringsAsFactor=F)
#targets$array <- as.factor(targets$array)
#targets=subset(targets, targets$Cohort == "QLD")

#rgSet=load("RGset.Robject");
#detP = detectionP(rgSet) ## calculate detection p-values

# Change directory to the analysis directory and load all previously computed Robject
setwd(directory) 
load("RGset.RObject");
rgSet=RGset
load("detP.RObject")

gmSetSwFlt = load("gmSetSwFlt.RObject") # Methylset values after SWAN normalization
gmSetQFlt = load("gmSetQFlt.RObject")# M values after SWAN normalization
mValsSw = load("mValsSw.RObject")# M values after SWAN normalization
mValsSq = load("mValsSq.RObject")# M values after SWAN normalization


#Density plot of values
#pd <- pData(RGset)

qcReport(RGset, pdf = "qcReportminfi.pdf")

jpeg("01_densityplot.jpg")
densityPlot(rgSet, main = "Beta", xlab = "Beta")
dev.off()


jpeg("02_densitybeanplot.jpg")
densityBeanPlot(rgSet)
dev.off()


jpeg("03_plotsex.jpg")
plotSex(rgSet)
dev.off()