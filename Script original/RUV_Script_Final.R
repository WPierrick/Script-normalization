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
source("~/ibscratch/MND/ExperimentalDesign/QLD/Script/utils.R")
source("~/ibscratch/MND/ExperimentalDesign/QLD/Script/functions.R")
ann.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
Xreact = read.csv(file="~/ibscratch/MND/ExperimentalDesign/Combined/IDAT_Chn_Qld/48639-non-specific-probes-Illumina450k.csv", 
                  stringsAsFactors=FALSE)

## Set directory
data="~/ibscratch/MND/ExperimentalDesign/Combined/IDAT_Chn_Qld"
directory="~/ibscratch/MND/ExperimentalDesign/Combined/RUV_Qld_V2"

#Set to data directory
setwd(data)

### Load Combined target data
sample_file="Chn_Qld_target.csv"

#twin data
twin="~/ibscratch/MND/ExperimentalDesign/QLD/20141003_Wray/Analysis2/DuplicateTwin.txt"
twin=read.table(twin, header=T)

## Specify covariates for which methylation data should be adjusted
covariates=c("predicted_sex","horvath_predicted_age","plate","array", 
             "array_pos","CD8T","CD4T","NK","Bcell","Mono","Gran")

## read in targets and idats
targets=read.csv(sample_file, header=T, stringsAsFactor=F)
targets$array <- as.factor(targets$array)
targets=subset(targets, targets$Cohort == "QLD")

rgSet=read.450k.exp(base=data, targets=targets,recursive=TRUE);
detP = detectionP(rgSet) ## calculate detection p-values

# Change directory to the analysis directory
setwd(directory) 
save(rgSet, file="rgSet.RObject")

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

## Remove poor quality probes
keepProbes = rowSums(detP < 0.01) == ncol(detP) 
mSetSwFlt = mSetSw[keepProbes,]
gmSetQFlt = gmSetQ[keepProbes,]
gmSetQFlt = gmSetQFlt[match(featureNames(mSetSwFlt),featureNames(gmSetQFlt)),]

## Remove probes with SNPs at CpG or single base extension (SBE) site
gmSetSwFlt = mapToGenome(mSetSwFlt)
gmSetSwFlt = dropLociWithSnps(gmSetSwFlt, snps = c("CpG", "SBE"))
gmSetQFlt = dropLociWithSnps(gmSetQFlt, snps = c("CpG", "SBE"))
gmSetQFlt = gmSetQFlt[match(featureNames(gmSetSwFlt),featureNames(gmSetQFlt)),]

## Remove cross-reactive probes
noXreact = !(featureNames(gmSetSwFlt) %in% Xreact$TargetID) 
gmSetSwFlt = gmSetSwFlt[noXreact,] 
noXreact = !(featureNames(gmSetQFlt) %in% Xreact$TargetID) 
gmSetQFlt = gmSetQFlt[noXreact,] 

## Remove sex shromosome probes
autosomes = !(featureNames(gmSetSwFlt) %in% ann.450k$Name[ann.450k$chr %in% c("chrX","chrY")])
gmSetSwFlt = gmSetSwFlt[autosomes,]
autosomes = !(featureNames(gmSetQFlt) %in% ann.450k$Name[ann.450k$chr %in% c("chrX","chrY")])
gmSetQFlt = gmSetQFlt[autosomes,]

## Get M-values
mValsSw = getM(gmSetSwFlt)
mValsSq = getM(gmSetQFlt)

save(gmSetSwFlt, file="gmSetSwFlt.RObject")
save(gmSetQFlt, file="gmSetQFlt.RObject")
save(mValsSw, file="mValsSw.RObject")
save(mValsSq, file="mValsSq.RObject")

##Remove individuals that does not have case and control information & one individual from a twin pair
targets <- targets[complete.cases(targets$CaseControl),]
targets <- targets[-which(targets$Basename %in% twin$Basename),]
targets$array <- as.factor(targets$array)
targets$Disease <- ifelse(targets$CaseControl=="Cases",1,0)
targets$Disease <- as.factor(targets$Disease)

##Remove these individuals from M-values
#mValsSw <- mValsSw[,which(colnames(mValsSw) %in% targets$Basename)]
mValsSq <- mValsSq[,which(colnames(mValsSq) %in% targets$Basename)]

#Make sure the IDs in targets and M values data in the same order
#all(targets$Basename == colnames(mValsSw))
all(targets$Basename == colnames(mValsSq))

#Adjust the probes for covariates
residuals=apply(mValsSq,1, function(x) batch_lm(as.numeric(x),targets[,covariates]))
residuals=t(residuals)
save(residuals, file="residuals.RObject")

#Inverse normal transformation
invN=apply(residuals,1, function(x) InvNorm(x))
invN=t(invN)
save(invN, file="invN.RObject")

#Change invN into mValsSq
mValsSq <- invN
colnames(mValsSq) <- targets$Basename
all(targets$Basename == colnames(mValsSq))

dmps = list(swan=list(),sqn=list())
Disease = factor(targets$Disease)
design = model.matrix(~Disease) ## Design matrix

#Compare different analysis model
##Linear model
assoc_results = apply(mValsSw,1,function(x) case_ctrl_assoc(x, targets$Disease));
assoc_results=as.data.frame(matrix(unlist(assoc_results), ncol=5, byrow=T));
colnames(assoc_results)=c("adjR2","effect","se","t" ,"pval");
assoc_results$Probe=rownames(mValsSw)
save(assoc_results, file="assoc_results.RObject");

p <- assoc_results$pval
png("LM.png")
qqplot(p,"LM")
dev.off()

##limma
dmps$sqn$limma = limmaFit(data=mValsSw,design=design,coef=2)
limma <- dmps$sqn$limma
save(limma, file="limma.RObject")

p <- limma$P.Value
png("limma.png")
qqplot(p,"LIMMA")
dev.off()


#RUV
#negM = getNegs(rgSet) 
#negSq = negM[,match(colnames(mValsSw),colnames(negM))]
#negSq = rbind(mValsSw,negSq)

#ctl = rownames(negSq) %in% rownames(negM)
#dmps$sqn$ruv1 = ruvFit(data=negSq, design=design, ctl=ctl, coef=2, method="inv")
#ruv1 <- dmps$sqn$ruv1
#save(ruv1, file="ruv1_inv.RObject")

ctl = rownames(mValsSw) %in% rownames(limma)[limma$adj.P.Val > 0.80]
dmps$sqn$ruv2 = ruvFit(data=mValsSw, design=design, ctl=ctl, coef=2)
ruv2 <- dmps$sqn$ruv2
save(ruv2, file="ruv2_inv_p080.RObject")

p <- ruv2$p.ebayes
png("ruv2_p80.png")
qqplot(p,"ruv2_p80")
dev.off()

ctl = rownames(mValsSw) %in% rownames(limma)[limma$adj.P.Val > 0.50]
dmps$sqn$ruv2 = ruvFit(data=mValsSw, design=design, ctl=ctl, coef=2)
ruv2 <- dmps$sqn$ruv2
save(ruv2, file="ruv2_inv_p050.RObject")

p <- ruv2$p.ebayes
png("ruv2_p50.png")
qqplot(p,"ruv2_p50")
dev.off()

ctl = rownames(mValsSw) %in% rownames(limma)[limma$adj.P.Val > 0.20]
dmps$sqn$ruv2 = ruvFit(data=mValsSw, design=design, ctl=ctl, coef=2)
ruv2 <- dmps$sqn$ruv2
save(ruv2, file="ruv2_inv_p020.RObject")

p <- ruv2$p.ebayes
png("ruv2_p20.png")
qqplot(p,"ruv2_p20")
dev.off()

## SVA
design0 = model.matrix(~1, data=targets)
dmps$sqn$sva = svaFit(data=mValsSw,design=design,design0=design0,coef=2)
sva <- dmps$sqn$sva
save(sva, file="sva.RObject")

p <- sva$P.Value
png("sva.png")
qqplot(p,"SVA")
dev.off()

## ISVA
dmps$sqn$isva = DoISVA(data.m=mValsSw, pheno.v=targets$Disease, cf.m=NULL)
isva <- dmps$sqn$isva
save(isva, file="isva.RObject")

p <- isva$spv
png("isva.png")
qqplot(p,"ISVA")
dev.off()