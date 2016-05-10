################################################################################
## Original file name: utils.R
## Additional functions used in analyses for "RUV for 450k methylation arrays"
## Written by: Jovana Maksimovic
## Last modified: 21/01/2015
################################################################################

## Identify probes that discriminate between blood cell sub types
## Code from minfi pickCompProbes function that identifies probes associated with cell type
getCellTypeProbes <- function(mSet, cellTypes = NULL, numProbes = 50) {
    
    splitit <- function(x) {
        split(seq(along=x), x)
    }
    
    p <- getBeta(mSet)
    pd <- as.data.frame(pData(mSet))
    if(is.null(cellTypes)) {
        cellTypes <- unique(pd$CellType)
    } else {
        if(!all(cellTypes %in% pd$CellType))
            stop("elements of argument 'cellTypes' is not part of 'mSet$CellType'")
        keep <- which(pd$CellType %in% cellTypes)
        pd <- pd[keep,]
        p <- p[,keep]
    }
    ## make cell type a factor
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    # get fstats
    ffComp <- rowFtests(p, pd$CellType)
    prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range")
    tIndexes <- splitit(pd$CellType)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })
    
    ## take numProbes up and numProbes down
    probeList <- lapply(tstatList, function(x) {
        y <- x[x[,"p.value"] < 1e-8,]
        yUp <- y[order(y[,"dm"], decreasing=TRUE),]
        yDown <- y[order(y[,"dm"], decreasing=FALSE),]
        c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
    })
    
    trainingProbes <- unlist(probeList)
    p <- p[trainingProbes,]
    
    return(p)
}

## Create ROC curve plot for list of p-values
ROCplot = function(rocList, truthVec, xlim=c(0,1), ylim=c(0,1), colours=rainbow(length(rocList)), main="ROC Plot", 
                   lwd=1, lty=1, legend=NULL, legendPos="bottomright", legendCex=1,legendLwd=lwd,legendLty=lty,
                   legendCols=colours,cex.lab=1,cex.axis=1,cex.main=1,maxP=NULL,xlab="False positive rate",
                   ylab="True positive rate",yaxt=c("default","none"),xaxt=c("default","none")){
    
    yaxt = match.arg(yaxt)
    xaxt = match.arg(xaxt)
    
    downsamp = NULL
    pred <- prediction(1-rocList[[1]]*10^100, truthVec)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
    
    fpr = unlist(perf@x.values)
    tpr = unlist(perf@y.values)
    #browser()
    if(!is.null(maxP)){
        downsamp = sort(sample(1:length(fpr), maxP))
        fpr = fpr[downsamp]
        tpr = tpr[downsamp]
    }
        
    plot(fpr, tpr, col="white", lwd=lwd[1], lty=lty[1], main=main, 
         xlim=xlim, ylim=ylim, axes=FALSE, xaxt="n", yaxt="n", ylab=ylab, 
         xlab=xlab,cex.lab=cex.lab,cex.main=cex.main)
    abline(0,1,lty=2)
    
    if(length(lty) == 1) lty=rep(lty,length(rocList))
    if(length(lwd) == 1) lwd=rep(lwd,length(rocList))
    
    for(i in 1:length(rocList)){
        
        if(!is.null(rocList[[i]])){
            pred <- prediction(1-rocList[[i]]*10^100, truthVec)
            perf <- performance(pred, measure = "tpr", x.measure = "fpr")
            
            fpr = unlist(perf@x.values)
            tpr = unlist(perf@y.values)
            
            if(!is.null(downsamp)){
                fpr = fpr[downsamp]
                tpr = tpr[downsamp]
            }
            
            par(new=TRUE)
            plot(fpr, tpr, type="l", col=colours[i], axes=FALSE, 
                 xaxt="n", yaxt="n", ylab="", xlab="", 
                 lwd=lwd[i], lty=lty[i], xlim=xlim, ylim=ylim)
        }
    }
    
    if(xaxt == "default") axis(1,at=seq(xlim[1],xlim[2],by=0.05),labels=seq(xlim[1],xlim[2],by=0.05), las=1, cex.axis=cex.axis)
    if(yaxt == "default") axis(2, at=seq(ylim[1],ylim[2],by=0.2), labels=seq(ylim[1],ylim[2],by=0.2), las=1, cex.axis=cex.axis)
    box()
    
    if(!is.null(legend)){
        legend(legendPos, legend=legend, cex=legendCex, col=legendCols, lty=legendLty, lwd=legendLwd,bg="white")
    }
    
    return()
}

## Prep vector of p-values for ROC analysis by ensuring order of p-values same as truth set
ROCprep <- function(pvals,truthSet){
    
    return(pvals[match(names(truthSet), names(pvals))])
}

## Run limma analysis and return top table
limmaFit = function(data,design,coef){
    
    fit = lmFit(data,design)
    fit2 = eBayes(fit)
    return(topTable(fit2,coef=coef,num=Inf))
    
}

## Run SVA analysis and return top table
svaFit = function(data,design,design0,coef){
    
    svobj = sva(data,design,design0)
    modSv = cbind(design,svobj$sv)
    fit = lmFit(data,modSv)
    fit2 = eBayes(fit)
    return(topTable(fit2,coef=coef,num=Inf))

}

## Run ComBat analysis and return top table
combatFit = function(data,design,batch,coef){
    
    combat_edata = ComBat(dat=data, batch=batch, mod=design, numCovs=NULL, par.prior=TRUE)
    
    fit = lmFit(combat_edata,design)
    fit2 = eBayes(fit)
    return(topTable(fit2,coef=coef,num=Inf))
}

## Run RUV analysis and return top table
ruvFit = function(data,design,ctl,coef,method=c("inv", "rinv", "ruv4", "ruv2")){
    
    method <- match.arg(method)
    fit = RUVfit(data=data, design=design, coef=coef, ctl=ctl, method=method, randomization=TRUE)
    fit2 = RUVadj(fit)
    return(topRUV(fit2, num=Inf))
    
}

## Get AUCs for list of p-values
getAUC <- function(rocList, truthSet){
    
    AUC = vector("numeric", length(rocList))
    
    for(i in 1:length(rocList)){
    
        pred <- prediction(1-rocList[[i]]*10^100, truthSet)
        perf <- performance(pred, measure = "auc")
        AUC[i] = unlist(perf@y.values)
    }
    
    return(AUC)
}

## Extract Illumina negative controls from minfi RGset object
## Return log2 ratio of green to red intensity as for Type II probes
getNegs <- function(rgSet){
    
    ctrls = getProbeInfo(rgSet, type = "Control")
    M.Neg = getGreen(rgSet)[ctrls$Address[ctrls$Type == "NEGATIVE"], ]
    U.Neg = getRed(rgSet)[ctrls$Address[ctrls$Type == "NEGATIVE"], ]
    
    M.Neg[M.Neg == 0] = min(M.Neg[M.Neg != 0])
    U.Neg[U.Neg == 0] = min(U.Neg[U.Neg != 0])
    
    log2(M.Neg/U.Neg)
}

## Subsample from large dataset by first sampling chips and then samples to achieve desired numbers
## Perform DM analysis using various methods on subsampled data
## Return list of top tables (length=reps) 
subSample <- function(reps,num,targets,raw,sw,fn,negM){
    
    chips = unique(targets$Chip)
    results = vector("list",reps)
    
    for(i in 1:reps){
        
        cat("\n****Round",i,"\n")
        
        set.seed(100+i)
        chip = sample(chips,1)
        tmp1 = targets[targets$Chip == chip, ]
        
        while(sum(table(tmp1$group) >= num) != 2){
            set.seed(100+i)
            chip = c(chip, sample(chips[!(chips %in% chip)],1))
            tmp1 = targets[targets$Chip %in% chip,]
        }
        
        set.seed(100+i)
        normal = as.character(sample(tmp1$Sample[tmp1$group == "Solid Tissue Normal"], num))
        cancer = as.character(tmp1$Sample[tmp1$group == "Primary solid Tumor" &
                                              tmp1$Chip %in% tmp1$Chip[tmp1$Sample %in% normal]])
        
        set.seed(100+i)
        if(length(cancer) < num){
            cancer = c(cancer, sample(as.character(tmp1$Sample[tmp1$group == "Primary solid Tumor" &
                                                                   !tmp1$Sample %in% cancer]), 
                                      (num-length(cancer))))
        } else {
            cancer = sample(cancer, num)
        }
        
        selected = targets$Sample %in% c(normal, cancer)

        tmp2 = list(swan=list(),fnorm=list())
        swSub = sw[,selected]
        fnSub = fn[,selected]
        
        status = factor(targets$group[selected])
        design = model.matrix(~status)
        
        cat("Limma\n")
        tmp2$swan$limma = limmaFit(data=swSub,design=design,coef=2)
        tmp2$fnorm$limma = limmaFit(data=fnSub,design=design,coef=2)
        cat("RUV\n")
        ctl = rownames(swSub) %in% rownames(tmp2$swan$limma)[ceiling(nrow(tmp2$swan$limma)*0.5):nrow(tmp2$swan$limma)]
        tmp2$swan$ruv = try(ruvFit(data=swSub, design=design, ctl=ctl, coef=2),TRUE)
        ctl = rownames(fnSub) %in% rownames(tmp2$fnorm$limma)[ceiling(nrow(tmp2$fnorm$limma)*0.5):nrow(tmp2$fnorm$limma)]
        tmp2$fnorm$ruv = try(ruvFit(data=fnSub, design=design, ctl=ctl, coef=2),TRUE)
        
        negSw = negM[,match(colnames(swSub),colnames(negM))]
        negSw = rbind(swSub,negSw)
        ctl = rownames(negSw) %in% rownames(negM)
        tmp2$swan$ruv1 = ruvFit(data=negSw, design=design, ctl=ctl, coef=2, method="inv")
        ctl = rownames(swSub) %in% rownames(tmp2$swan$ruv1)[dmps$swan$ruv1$p.ebayes.BH > 0.2]
        tmp2$swan$ruv2 = ruvFit(data=swSub, design=design, ctl=ctl, coef=2)
        
        negFn = negM[,match(colnames(fnSub),colnames(negM))]
        negFn = rbind(fnSub,negFn)
        ctl = rownames(negFn) %in% rownames(negM)
        tmp2$fnorm$ruv1 = ruvFit(data=negFn, design=design, ctl=ctl, coef=2, method="inv")
        ctl = rownames(fnSub) %in% rownames(tmp2$fnorm$ruv1)[tmp2$fnorm$ruv1$p.ebayes.BH > 0.2]
        tmp2$fnorm$ruv2 = ruvFit(data=fnSub, design=design, ctl=ctl, coef=2)
        
        cat("ISVA\n")
        tmp2$swan$isva = try(DoISVA(data.m=swSub, pheno.v=targets$group[selected], cf.m=NULL),TRUE)
        tmp2$fnorm$isva = try(DoISVA(data.m=fnSub, pheno.v=targets$group[selected], cf.m=NULL),TRUE)
        cat("SVA\n")
        design0 = model.matrix(~1, data=targets[selected,])
        tmp2$swan$sva = try(svaFit(data=swSub,design=design,design0=design0,coef=2),TRUE)
        tmp2$fnorm$sva = try(svaFit(data=fnSub,design=design,design0=design0,coef=2),TRUE)
        
        results[[i]] = tmp2
    }
    
    return(results)
}

plotAUCs <- function(auc,max,reps,colours,ylim,key,scr,yaxt,ylab,samps,showKey){
    
    screen(scr)
    #mar=c(bottom, left, top, right)
    par(mar=c(5,0,5,0),oma=c(0,5,0,3),xpd=NA)
    
    for(i in 1:(length(auc)/reps)){
        
        if(i == 1){
            
            plot(1:reps,auc[1:reps],col=colours[i],ylim=ylim,pch=ifelse(auc[1:reps] %in% max,19,1),
                 ylab=ylab,yaxt=yaxt,xaxt="n",xlab="")
            axis(1,at=1:reps)
            if(showKey) legend("bottomright",legend=key,col=colours,pch=1)
            title(xlab="Sampling round")
            title(paste(samps, "samples per group"),line=1)
        } else {
            
            sapply(1:3, function(j) {
                points(1:reps,auc[((reps*j)+1):(reps*(j+1))],col=colours[j+1],
                       pch=ifelse(auc[((reps*j)+1):(reps*(j+1))] %in% max,19,1),yaxt=yaxt)
            })
        }
    }
}


