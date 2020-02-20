assess <- function(sim.out, alpha.type=c("fdr","pval"), alpha.level=0.05,
                   stratify.type=c("abundance", "prevalence"),stratify,stratify.num=10) {
  
  alpha.type = match.arg(alpha.type)
  stratify.type = match.arg(stratify.type)
  
  
  if(stratify.type=='prevalence'){
    if(missing(stratify)) stratify = c(0,0.2,0.4,0.8,1)
  }else if(stratify.type=='abundance'){
    maxbase2=round(max(log(sim.out$x.abunds)),2)
    if(missing(stratify)) stratify = c(0,2^(1:maxbase2))
  }
  
  n1s = sim.out$n1s
  n2s = sim.out$n2s
  nTaxa = sim.out$params$nTaxa
  params = sim.out$params
  DAids = sim.out$DAids
  lfcs.mu = sim.out$lfcs.mu
  nsims = dim(sim.out$pvalue)[3]
  
  stratify.name=paste("(",stratify[-length(stratify)],',',stratify[-1],"]",sep='')
  TP.stratify = FP.stratify = FDR.stratify = TPR.stratify = FPR.stratify = taxon.stratify=array(NA,dim=c(length(stratify) - 1,length(n1s), nsims),dimnames = list(stratify.name=stratify.name,samplesize=paste(n1s,n2s,sep = '/'),sim=1:nsims))
  TP = FP = FDR = TPR =FPR = array(NA,dim=c(length(n1s), nsims),dimnames=list(samplesize=paste(n1s,n2s,sep = '/'),sim=1:nsims))
  
  
  for(i in 1:nsims) {
    for(j in 1:length(n1s)) {
      
      nn1 = n1s[j]
      nn2 = n2s[j]
      DAid = sim.out$DAids[[i]]
      lfc = sim.out$lfcs.mu[[i]]
      Zg= rep(0, nTaxa)
      Zg[DAid] = 1
      
      if(stratify.type == "abundance") {
        x = sim.out$x.abunds[,j,i]
      }else  if(stratify.type == "prevalence") {
        x = sim.out$x.prevs[,j,i]
      }
      xgr = cut(x, stratify)
      stratify.keep=which(table(xgr)>=stratify.num)
      ix.keep = seq(length(x))[xgr %in% levels(xgr)[stratify.keep]]
      xgr = cut(x[ix.keep], stratify[c(stratify.keep,max(stratify.keep)+1)])
      ix.keep2=match(names(table(xgr)),stratify.name)
      taxon.stratify[ix.keep2,j,i]=table(xgr)
      Zg = Zg[ix.keep]
      
      if(alpha.type == "pval"){
        x = sim.out$pvalue[ix.keep,j,i]
      }else {
        #x = sim.out$fdr[ix.keep,j,i]
        pval = sim.out$pvalue[ix.keep,j,i]
        x = p.adjust(pval, "BH")
      }
      
      assess0.out = assess0(x, alpha.level, Zg, xgr=xgr)
      
      TP[j,i]=assess0.out$TP
      TP.stratify[ix.keep2,j,i] = assess0.out$TP.stratify
      TPR[j,i]=assess0.out$TPR
      TPR.stratify[ix.keep2,j,i] = assess0.out$TPR.stratify
      
      FP[j,i]=assess0.out$FP
      FP.stratify[ix.keep2,j,i] = assess0.out$FP.stratify
      FPR[j,i]=assess0.out$FPR
      FPR.stratify[ix.keep2,j,i] = assess0.out$FPR.stratify
      
      FDR[j,i] = assess0.out$FDR
      FDR.stratify[ix.keep2,j,i] = assess0.out$FDR.stratify
      
    }
  }
  
  output <- list(TP=TP,TP.stratify=TP.stratify, 
                 TPR=TPR,TPR.stratify=TPR.stratify, 
                 FP=FP,FP.stratify=FP.stratify, 
                 FPR=FPR,FPR.stratify=FPR.stratify, 
                 FDR=FDR,FDR.stratify=FDR.stratify, 
                 taxon.stratify=taxon.stratify,
                 alpha.type=alpha.type, alpha.level=alpha.level,
                 stratify.type=stratify.type, stratify.name=stratify.name,
                 n1s=n1s,n2s=n2s)
  
}




assess0 <- function(p, alpha.level, Zg, xgr){
  
  ix.P = p <= alpha.level
  P = sum(ix.P) # P=TP+FP
  P.stratify = tapply(ix.P, xgr, sum)
  
  #TP, TPR (power)
  id.TP = Zg==1
  TP.stratify = tapply(p[id.TP] <= alpha.level, xgr[id.TP], sum)
  TP.stratify[is.na(TP.stratify)] = 0
  TP=sum(TP.stratify,na.rm=TRUE)
  TPR=TP/sum(id.TP)
  TPR.stratify=as.vector(TP.stratify/table(xgr[id.TP]))
  
  #FP, FPR (type I error)
  id.TN = Zg==0
  FP.stratify = tapply(p[id.TN] <= alpha.level, xgr[id.TN], sum)
  FP.stratify[is.na(FP.stratify)] = 0
  FP=sum(FP.stratify,na.rm=TRUE)
  FPR=FP/sum(id.TN)
  FPR.stratify=as.vector(FP.stratify/table(xgr[id.TN]))
  
  #FDR
  P.stratify=apply(P.stratify,1,function(x) max(1,x))
  FDR.stratify = FP.stratify / P.stratify
  #FDR = sum(FP.stratify, na.rm=TRUE) / max(1,P)
  FDR=FP/max(1,P)
  
  list(alpha.level=alpha.level,
       TP=TP,TP.stratify=TP.stratify,
       TPR=TPR,TPR.stratify=TPR.stratify,
       FP=FP,FP.stratify=FP.stratify,
       FPR=FPR,FPR.stratify=FPR.stratify,
       FDR=FDR,FDR.stratify=FDR.stratify)
}



summaryAssess<- function(assess.out,assess.type=c('overall','stratify'),is.print=TRUE) {
  
  assess.type=match.arg(assess.type)
  
  if(assess.type=='overall'){
    
    TP.m=rowMeans(assess.out$TP,na.rm=T)
    TPR.m=rowMeans(assess.out$TPR,na.rm=T)
    FP.m=rowMeans(assess.out$FP,na.rm=T)
    FPR.m=rowMeans(assess.out$FPR,na.rm=T)
    FDR.m= rowMeans(assess.out$FDR,na.rm=T)
    
    res.m=cbind(assess.out$n1s, assess.out$n2s,assess.out$alpha.level, FDR.m, TP.m, FP.m, TPR.m, FPR.m)
    if(assess.out$alpha.type == "fdr") {
      colnames(res.m)= c("G1","G2","FDR.nominal","FDR", "TP", "FP", "TPR","FPR")
    }else{
      colnames(res.m)= c("G1","G2","Type1Error.nominal","FDR", "TP", "FP", "TPR","FPR")
    }
    
  }else if(assess.type=='stratify'){
    
    res.m=res.s=array(NA,dim=c(length(assess.out$stratify.name),length(assess.out$n1s), 6),
                      dimnames = list(stratify.name=assess.out$stratify.name,
                                      samplesize=paste(assess.out$n1s,assess.out$n2s,sep = '/'),
                                      crit=c('TP','TPR','FP','FPR','FDR','Taxon')))
    
    TP.m=apply(assess.out$TP.stratify,c(1,2), mean, na.rm=TRUE)
    TPR.m=apply(assess.out$TPR.stratify,c(1,2), mean, na.rm=TRUE)
    FP.m=apply(assess.out$FP.stratify,c(1,2), mean, na.rm=TRUE)
    FPR.m=apply(assess.out$FPR.stratify,c(1,2), mean, na.rm=TRUE)
    FDR.m=apply(assess.out$FDR.stratify,c(1,2), mean, na.rm=TRUE)
    taxon.m=apply(assess.out$taxon.stratify,c(1,2), mean, na.rm=TRUE)
    res.m[,,1]=TP.m;res.m[,,2]=TPR.m; res.m[,,3]=FP.m;res.m[,,4]=FPR.m; res.m[,,5]=FDR.m; res.m[,,6]=taxon.m
    
    
    TP.s=apply(assess.out$TP.stratify,c(1,2), sd, na.rm=TRUE)
    TPR.s=apply(assess.out$TPR.stratify,c(1,2), sd, na.rm=TRUE)
    FP.s=apply(assess.out$FP.stratify,c(1,2), sd, na.rm=TRUE)
    FPR.s=apply(assess.out$FPR.stratify,c(1,2), sd, na.rm=TRUE)
    FDR.s=apply(assess.out$FDR.stratify,c(1,2), sd, na.rm=TRUE)
    taxon.s=apply(assess.out$taxon.stratify,c(1,2), sd, na.rm=TRUE)
    res.s[,,1]=TP.s;res.s[,,2]=TPR.s; res.s[,,3]=FP.s;res.s[,,4]=FPR.s; res.s[,,5]=FDR.s; res.s[,,6]=taxon.s
    
  }
  
  if(is.print) print(signif(res.m,3))
  
  if(assess.type=='overall'){
    return((signif(res.m,3)))
  }else  if(assess.type=='stratify'){
    return(list(m=signif(res.m,3),s=signif(res.s,3)))
  }
  
}



