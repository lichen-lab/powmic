
.metagenomeSeq<-function(dat) {
  mgs=newMRexperiment(counts=dat$counts)
  mgs=metagenomeSeq::cumNorm(mgs)
  mod=model.matrix(~dat$designs)
  settings=zigControl(verbose = F)
  fit=fitZig(mgs,mod,control=settings)
  res = MRfulltable(fit,number = nrow(assayData(mgs)$counts))
  ix=as.numeric(rownames(res))
  pval=res$pvalues
  pval[is.na(pval)]=1
  fdr=res$adjPvalues
  fdr[is.na(fdr)]=1
  res=data.frame(ix=ix,pval=pval,fdr=fdr)
  res=res[order(res$ix,decreasing = F),]
  res
}




.wilcox <- function(dat) {
  pval=apply(dat$counts,1,function(x) {wilcox.test(x[dat$designs==0],x[dat$designs==1])$p.value})
  pval[is.na(pval)]=1
  fdr=p.adjust(pval,method="fdr")
  res=data.frame(ix=1:length(pval),pval=pval,fdr=fdr)
  res
}



.edgeR <- function(dat) {
  d=DGEList(counts=dat$counts, group=dat$designs)
  d=edgeR::calcNormFactors(d)
  d=estimateCommonDisp(d)
  d=estimateTagwiseDisp(d)
  fit=exactTest(d)
  res=topTags(fit, n=nrow(dat$counts))
  res=as.data.frame(res)
  ix=as.numeric(rownames(res))
  pval=res[,'PValue']
  fdr=res[,'FDR']
  res=data.frame(ix=ix,pval=pval,fdr=fdr)
  res=res[order(res$ix),]
  res
}


.DESeq2 <- function(dat) {
  cond=factor(dat$designs)
  d=DESeqDataSetFromMatrix(dat$counts, DataFrame(cond), ~ cond)
  sf.tss=colMeans(dat$counts);sf.tss=sf.tss/median(sf.tss)
  sizeFactors(d)=sf.tss
  d=DESeq(d, quiet=TRUE)
  res=results(d)
  pval=res$pvalue
  pval[is.na(pval)]=1
  fdr=res$padj
  fdr[is.na(fdr)]=1
  res=data.frame(ix=1:length(pval), pval=pval, fdr=fdr)
  res
}


