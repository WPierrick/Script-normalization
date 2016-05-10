
# 2014-09-29
# Sonia Shah
# s.shah1@uq.edu.au



## logit transform
logit_transform=function(b) {
    b=log(b/(1-b))
    return(as.data.frame(b))
}


## Adjust for batch effect
batch_lm=function(b, mybatch) {
    reg=lm(b~., data=mybatch, na.action=na.exclude)
    residuals=residuals(reg, type="working")
    return(residuals)
}

## Inverse normal transformation
## Function written by Zhihong Zhu
InvNorm=function(x) {
    return(qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x))))
}

## Zscore for phenotype after adjusting for covariates
## Function written by Zhihong Zhu
ZScore=function(y,x) {
    sd=sd(y, na.rm=T)
    mean=mean(y,na.rm=T)
    b=coefficients(summary(lm(y~., data=x)))
    e=(y-as.matrix(x)%*%b[-1,])[,1]
    return((e-mean(e,na.rm=T))/sd(e,na.rm=T))
}


## association analysis function
## No covariates as both trait and methylation data should have been adjusted for covariates separately
assoc_analysis=function(b, t) {
  
    reg=summary(lm(t ~ b))
    effect=reg$coefficients[2,"Estimate"]
    se=reg$coefficients[2,"Std. Error"]
    pval=reg$coefficients[2,"Pr(>|t|)"]
    return(list(effect=effect, se=se, pval=pval))
}


## case ctrl analysis
#assoc function
case_ctrl_assoc=function(b, status) {
  reg=summary(lm(b~status, na.action=na.exclude))
  effect=reg$coefficients[2,"Estimate"]
  r2=reg$adj.r.squared
  se=reg$coefficients[2,"Std. Error"]
  t=reg$coefficients[2,"t value"]
  pval=reg$coefficients[2,"Pr(>|t|)"]
  
  return(list(adjR2=r2, effect=effect, se=se, t=t, pval=pval))
  
}



## Run EWAS
run_ewas= function(b, pheno, traits, covars) {
    
    for (t in traits) {
        
        system(paste("mkdir ./Results/", t, sep=""))
        
        ## Run EWAS
        ## Generate sex and age-adjusted zscore values for phenotype 
        zscore=ZScore(pheno[,t], pheno[,covars])
        
        ## No covariates passed to assoc_analysis as both trait and methylation data should have been adjusted for covariates separately
        assoc_results=apply(b,1, function(x) assoc_analysis(x, zscore))
        assoc_results=data.frame(matrix(unlist(assoc_results), nrow=nrow(b), byrow=T))
        colnames(assoc_results)=c("effect","SE","P")
        assoc_results$Probe=rownames(b)
        assoc_results=assoc_results[,c("Probe","effect","SE","P")]
        write.table(assoc_results, file=paste("./Results/", t , "/" , t ,"_EWAS_results.txt", sep=""), sep="\t", quote=F, row.names=F)
        
        jpeg(paste("./Results/", t ,"/qqplot_", t , ".jpg",sep=""))
        qq(assoc_results$P)
        dev.off()
    }
}





## Get uncorrelated probes
## dat should be a data.frame
prune_probes = function(results, dat, pval, cor_threshold, window) {
    
    library(reshape)
    
    # Read in probe annotation
    #annot=read.csv("./Annovar_hg19_probe_annot.csv", stringsAsFactor=F, header=T)
    #load("Illumina450k_probe_annot.RObject")
    #results=merge(results, probe_annot, by.x="Probe",by.y="probe_id", all.x=T, all.y=F)
    #rownames(results)=results$Probe
    
    # sort ewas results by chr:pos
    results=sort_df(results, vars=c("chromosome","Coordinate_37"))
    
    # Get probes that pass the specfied p-value threshold
    results=results[which(results$P <= pval),]
    pruned_results=results
    
    # Sliding window approach
    x=1
    start_chr=results[x,"chromosome"]
    start_pos=results[x,"Coordinate_37"]
    end_pos=start_pos + window
    corprobes=c()
    
    probes_to_remove=c()
    dat=dat[results$Probe,]
    temp=dat
    
    while(x <= nrow(dat)) {
            
        # get all probes in 250bp window
        probes=results$Probe[which(results$chromosome==start_chr & results$Coordinate_37>=start_pos & results$Coordinate_37<=end_pos)]
        
        if (length(probes)>1) {
        
        
            temp=dat[probes,]
   
         
             # pairwise correlation with most significant probe in window
            temp_signif=results[probes,]
            temp_signif=sort_df(temp_signif, vars="P")
            top_probe=temp_signif$Probe[1]
            
            
            # if there are more than 2 probes in addition to the top probe in the window
            if (length(probes)>2) {
                correl=abs(cor(x=t(temp[top_probe,]), y=t(temp[-which(rownames(temp)==top_probe),]), use="pairwise.complete.obs"))
                corprobes=colnames(correl)[which(correl>cor_threshold)]
        
            # if there is only one probe in addition to the top probe in the window
            } else {
                 correl=abs(cor(x=as.double(temp[top_probe,]), y=as.double(temp[-which(rownames(temp)==top_probe),]), use="pairwise.complete.obs"))
                if (correl>cor_threshold) { corprobes = rownames(temp)[which(rownames(temp)!=top_probe)]}
                
            }
                # remove probes with correlation greater than specified threshold
          

            
            if (length(corprobes)!=0) {
                probes_to_remove=c(probes_to_remove, corprobes)
                
                p=which(rownames(dat) %in% probes_to_remove)
                
                if (length(p) >0) {
                    dat=dat[-p,]
                    results=results[-p,]
                }
            }
            

        }
        
        x=x+1
        start_chr=results[x,"chromosome"]
        start_pos=results[x,"Coordinate_37"]
        end_pos=start_pos + window
    } # end while loop
    
    if (length(probes_to_remove)!=0) {
        pruned_results=pruned_results[-which(pruned_results$Probe %in% probes_to_remove),]
    } 
    
    return(pruned_results)
}




## Generate Height EWAS predictor
ewas_score= function(d,  effects) {
   
        score=c()
        d=d[effects$Probe,]
        x=rowSums(t(effects$effect*d), na.rm=T)
    
        score$sampleID=names(x)
        score$ewas_score=x
        score=as.data.frame(score)
        score$sampleID=as.character(score$sampleID)
        return(score)
   
}



## Generate Height EWAS predictor
predict_trait= function(t, covars, mypheno, score, outfile) {
    
    sink(paste("./Results/",t,"/",outfile_prefix,"_",t,"_ewas_prediction_results.txt", sep=""))
    zscore=ZScore(mypheno[,t], mypheno[,covars])
        
    summary(lm(zscore ~ score))
    sink()

}

## QQ plot with lambda

qqplot <- function(x, title) {
      p1 <- x
      p2 <- sort(t(p1))
      n  <- length(p2)
      k  <- c(1:n)
      alpha   <- 0.05
      lower   <- qbeta(alpha/2,k,n+1-k)
      upper   <- qbeta((1-alpha/2),k,n+1-k)
      expect  <- (k-0.5)/n
      biggest <- ceiling(max(-log10(p2),-log10(expect)))
      lambd = median(qchisq(1-p1,1,lower.tail = TRUE)/qchisq(0.5, df=1,lower.tail = TRUE))
      shade <- function(x1, y1, x2, y2, color = col.shade) {
            n <- length(x2)
            polygon(c(x1, x2[n:1]), c(y1, y2[n:1]), border = NA, col = color)
      }
      plot(-log10(expect),-log10(p2),main=title,xlim=c(0,biggest),ylim=c(0,biggest),
           ylab=expression(paste("Observed  ","-",log[10],"(P)")), 
           xlab=expression(paste("Expected  ","-",log[10],"(P)")), type = "n")
      text(-0.15, y=biggest-0.3, paste("lambda = ", signif(lambd,3),sep=""), pos = 4, font=2)
      shade(-log10(expect),-log10(lower),-log10(expect),-log10(upper), color = "gray")
      abline(0,1,col="white",lwd=2)
      points(-log10(expect),-log10(p2), pch=20)
}



