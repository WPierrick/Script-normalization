### RUV Analysis for Parkinson disease methylation samples ###
### Modified by Pierrick Wainschtein from Beben Benyamin from Maksimovic et al 2015 ####

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
source("/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/utils.R")
source("/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/functions.R")
ann.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
Xreact = read.csv(file="/clusterdata/uqpwains/ibscratch/workdir/scripts/QCscript/48639-non-specific-probes-Illumina450k.csv", 
                  stringsAsFactors=FALSE)

## Set directory
data="/clusterdata/uqpwains/ibscratch/RAW-IDATS"
directory="/clusterdata/uqpwains/ibscratch/QC-IDATS"

#Set to data directory
setwd(data)

### Load Combined target data
sample_file="PD_Targets_09May2016.csv"

#twin data
#twin="~/ibscratch/MND/ExperimentalDesign/QLD/20141003_Wray/Analysis2/DuplicateTwin.txt"
#twin=read.table(twin, header=T)

## Specify covariates for which methylation data should be adjusted
covariates=c("SAMPLE.SET", "Gender", "Phenotype", "Sample.Well.for.Sample.Sheet", "Sentrix.Barcode", "Sample.Section")

## read in targets and idats
targets=read.csv(sample_file, header=T, stringsAsFactor=F)
targets$"X2D.MATRIX.ID" <- as.factor(targets$"X2D.MATRIX.ID")
#targets=subset(targets, targets$Cohort == "QLD")

#rgSet=read.450k.exp(base=data, targets=targets,recursive=TRUE);
load("RGset.Robject");
rgSet=RGset
detP = detectionP(rgSet) ## calculate detection p-values

# Change directory to the analysis directory
setwd(directory) 
save(detP, file="detP_raw.RObject")

#######################################################################################

## Examine mean detection p-values for all samples
pdf("01-detection-p-values.pdf", width=14)
barplot(apply(detP,2,mean), main="mean detection p-values", col=as.factor(targets$Source), xaxt="none")
abline(h=0.001,col="red")
dev.off()

#MDS plot on raw beta
png("01_MDS_plot.png")
mdsPlot(RGset, sampNames = targets$SAMPLE.SET, sampGroups = targets$SAMPLE.SET, 
        pch = 1, pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
        legendNCol, main = NULL)



png("01_bisMDS_plot.png")
par(mfrow=c(1,1))
plotMDS(RGset, labels=targets$Phenotype, col=as.integer(factor(targets$Phenotype)))
legend("topleft", legend=c("PD", "Control"), pch=16,cex=1.2,col=1:2)
dev.off()


png("01_bis2MDS_plot.png")
par(mfrow=c(1,1))
plotMDS(RGset)
dev.off()

# Density plot 
jpeg("01_densityplotgender.jpg")
densityPlot(RGset, sampGroups = targets$Gender, main = "Beta sorted by gender", xlab = "Beta")
dev.off()
jpeg("01_densityplotsampleset.jpg")
densityPlot(RGset, sampGroups = targets$SAMPLE.SET, main = "Beta by sample set", xlab = "Beta")
dev.off()
jpeg("01_densityplotphenotype.jpg")
densityPlot(RGset, sampGroups = targets$Phenotype, main = "Beta sorted by phenotype", xlab = "Beta")
dev.off()


#Get a Mset raw
MSet.raw <- preprocessRaw(RGset)
save(MSet.raw , file="MSet.raw.RObject")

## Pre-process the data after excluding poor quality samples
set.seed(100)
mSetSw = preprocessSWAN(rgSet[,apply(detP,2,mean) < 0.001]) ## probe-type normalisation only
gmSetQ = preprocessQuantile(rgSet[,apply(detP,2,mean) < 0.001]) ## probe-type and btw array normalisation
gmSetQ = gmSetQ[match(featureNames(mSetSw),featureNames(gmSetQ)),] ## ensure probes are ordered the same

## Remove poor quality samples from targets info and detection p-values  
targets = targets[apply(detP,2,mean) < 0.001,]
detP = detP[,apply(detP,2,mean) < 0.001]

save(targets, file="targets.RObject")
save(detP, file="detP.RObject")

pdf("02-detection-p-values-after-removal.pdf", width=14)
barplot(apply(detP,2,mean), main="Mean detection p-values after renoval", col=as.factor(targets$X2D.MATRIX.ID), xaxt="none")
abline(h=0.001,col="red")
dev.off()

#Identify the removed individuals
out <- which(!colnames(detP_raw)%in%colnames(detP))
IDremoved <- colnames(detP_raw)[out]
save(IDremoved, file="removedindiv.RObject")

## Remove poor quality probes
keepProbes = rowSums(detP < 0.001) == ncol(detP) 
mSetSwFlt = mSetSw[keepProbes,]
gmSetQFlt = gmSetQ[keepProbes,]
gmSetQFlt = gmSetQFlt[match(featureNames(mSetSwFlt),featureNames(gmSetQFlt)),]


# Density plot 
#jpeg("02_densityplotgender.jpg")
#densityPlot(detP, sampGroups = targets$Gender, main = "Beta sorted by gender", xlab = "Beta")
#dev.off()
#jpeg("02_densityplotsampleset.jpg")
#densityPlot(detP, sampGroups = targets$SAMPLE.SET, main = "Beta by sample set", xlab = "Beta")
#dev.off()
#jpeg("02_densityplotphenotype.jpg")
#densityPlot(detP, sampGroups = targets$Phenotype, main = "Beta sorted by phenotype", xlab = "Beta")
#dev.off()



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

save(gmSetSwFlt, file="gmSetSwFlt.RObject") # MethylSet with SWAN normalization
save(gmSetQFlt, file="gmSetQFlt.RObject") # MethylSet with quantile normalization
save(mValsSw, file="mValsSw.RObject") # M values with SWAN normalization
save(mValsSq, file="mValsSq.RObject") # M values  with quantile normalization


#Clustering
revmval  <- t(mValsSw)
cluster1 = hclust(dist(revmval, method = "euclidean"), method = "ward.D2")
png("03_cluster.png")
plot(cluster1,"Clustering")
dev.off()

# Effect of normalizing using SWAN
png("03_Normalization_plot.png")
par(mfrow=c(1,2))
plotBetasByType(MSet.raw [,1], main = "Raw")
plotBetasByType(gmSetSwFlt[,1], main = "SWAN")
dev.off()

##Remove individuals that does not have case and control information & one individual from a twin pair
targets <- targets[complete.cases(targets$Phenotype),]
targets$X2D.MATRIX.ID <- as.factor(targets$X2D.MATRIX.ID)
targets$Phenotype <- ifelse(targets$Phenotype=="PD",1,0)
targets$Phenotype <- as.factor(targets$Phenotype)

##Remove these individuals from M-values
#mValsSw <- mValsSw[,which(colnames(mValsSw) %in% targets$X2D.MATRIX.ID)]
#mValsSq <- mValsSq[,which(colnames(mValsSq) %in% targets$X2D.MATRIX.ID)]

#MDS plot on clean individuals
png("04_MDS_plot.png")
par(mfrow=c(1,1))
plotMDS(Mval, labels=targets$Phenotype, col=as.integer(factor(targets$Phenotype)))
legend("topleft", legend=c("PD", "Control"), pch=16,cex=1.2,col=1:2)
dev.off()

#Make sure the IDs in targets and M values data in the same order
all(targets$X2D.MATRIX.ID == colnames(mValsSw))
all(targets$X2D.MATRIX.ID == colnames(mValsSq))

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
colnames(mValsSq) <- targets$X2D.MATRIX.ID
all(targets$X2D.MATRIX.ID == colnames(mValsSq))

dmps = list(swan=list(),sqn=list())
Disease = factor(targets$Phenotype)
design = model.matrix(~Disease) ## Design matrix


#Compare different analysis model
##Linear model
assoc_results = apply(mValsSw,1,function(x) case_ctrl_assoc(x, targets$Phenotype));
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