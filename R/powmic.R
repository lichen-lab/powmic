
powmic <- function(n1s=c(20,40,60), n2s=c(20,40,60), nsims=10, params,
                    distrib=c('NB','ZINB','ZIP','DM'),
                    DAmethod=c("edgeR", "DESeq2","metaGenomeSeq","wilcox"), 
                    prevalence=0,verbose=TRUE) {
  
  distrib=match.arg(distrib)
  DAmethod = match.arg(DAmethod)
  
  setParams=switch(distrib,ZIP=setParams.ZIP, NB=setParams.NB, ZINB=setParams.ZINB, DM=setParams.DM)
  simCounts=switch(distrib,ZIP=simCounts.ZIP, NB=simCounts.NB, ZINB=simCounts.ZINB, DM=simCounts.DM)
  
  
  n1 = max(n1s)
  n2 = max(n2s)
  
  set.seed(params$sim.seed)
  pvalues = fdrs = x.prevs =x.abunds = array(NA,dim=c(params$nTaxa,length(n1s), nsims))
  DAids  = NULL
  
  for(i in 1:nsims) {
    if(verbose)
      cat("Simulation number", i, "\n")
    
    if(distrib=='NB'){
      params=setParams(nTaxa = params$nTaxa,p.DA=params$p.DA,Sigma=params$Sigma,
                           lmu0=params$lmu0,lphi=params$lphi0,
                           lfc.mu=params$lfc.mu,
                           depth=params$depth,
                           sim.seed=params$sim.seed+1)
      
    }else if(distrib=='ZINB'){
      params=setParams(nTaxa = params$nTaxa,p.DA=params$p.DA,Sigma=params$Sigma,
                           lmu0=params$lmu0,lphi=params$lphi0,lp0=params$lp0,
                           lfc.mu=params$lfc.mu,
                           depth=params$depth,
                           sim.seed=params$sim.seed+1)
    }else if(distrib=='ZIP'){
      params=setParams(nTaxa = params$nTaxa,p.DA=params$p.DA,Sigma=params$Sigma,
                           lmu0=params$lmu0,lp0=params$lp0,
                           lfc.mu=params$lfc.mu,
                           depth=params$depth,
                           sim.seed=params$sim.seed+1)
    }else if(distrib=='DM'){
      params=setParams(nTaxa = params$nTaxa,p.DA=params$p.DA,
                           lmu0=params$lmu0,lphi0=params$lphi0,
                           lfc.mu=params$lfc.mu,
                           depth=params$depth,libsize=params$libsize,
                           sim.seed=params$sim.seed+1)
    }
    
    sim.counts = simCounts(params, n1, n2)
    DAids[[i]] = sim.counts$DAid
    
    for(j in 1:length(n1s)) {
      
      nn1 = n1s[j]
      nn2 = n2s[j]
      
      idx = c(1:nn1, n1+(1:nn2))
      designs = sim.counts$designs[idx]
      x = sim.counts$counts[,idx]
      
      ss = rowMeans(x!=0)
      ix= ss>=prevalence
      x = x[ix,, drop=FALSE]
      
      dat=list(counts=x, designs=designs)
      if (DAmethod=="edgeR"){
        DA.res= suppressWarnings(.edgeR(dat))
      }else if (DAmethod=="DESeq2"){
        DA.res= suppressWarnings(.DESeq2(dat))
      }else if(DAmethod=="metaGenomeSeq"){
        DA.res= suppressWarnings(.metaGenomeSeq(dat))
      }else if (DAmethod=="wilcox"){
        DA.res= suppressWarnings(.wilcox(dat))
      }
      
      pval = fdr = numeric(nrow(x))+1
      x.abund=x.prev = numeric(nrow(x))
      pval[ix] = DA.res[, "pval"]
      fdr[ix] = DA.res[, "fdr"]
      sfs = colSums(dat$counts); sfs = sfs/median(sfs)
      x.abund[ix] = rowMeans(sweep(dat$counts,2,sfs,FUN="/"))
      x.prev[ix] = rowMeans(dat$counts!=0)
      
      pvalues[,j,i] = pval
      fdrs[,j,i] = fdr
      x.abunds[,j,i]=x.abund
      x.prevs[,j,i]=x.prev
      
    }
  }
  
  obj=list(pvalues=pvalues, fdrs=fdrs, x.prevs=x.prevs, x.abunds=x.abunds, DAids=DAids, n1s=n1s, n2s=n2s, params=params)
  structure(obj,class="powmic")
  
}





