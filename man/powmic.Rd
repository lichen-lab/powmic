\name{powmic}
\alias{powmic}
\title{Power assessment in microbiome case-control studies}
\description{
    Perform simulation-based power assessment for differential abundance analysis of metagenomic sequencing data. Microbiome data can be generated from distributions included Negative Binomial (NB), Zero-Inflated Negative Binomial (ZINB), Zero-Inflated Poisson Distribution (ZIP) and Dirichlet-Multinomial Distribution (DM). The parameters of these distributions as well as the microbial correlation network can be estimated from real datasets. The estimated parameters are used toe generate synthetic datasets for the power assessment. TPR and FDR are reported for different sample sizes. The users can choose a sample size with desired TPR and FDR in the study design.
}

\usage{
powmic(n1s=c(20,40,60), n2s=c(20,40,60), nsims=10, params,
                    distrib=c('NB','ZINB','ZIP','DM'),
                    DAmethod=c("edgeR", "DESeq2","metaGenomeSeq","wilcox"), 
                    prevalence=0,verbose=TRUE) 
}

\arguments{
  \item{n1s}{Sample sizes in case group}
  \item{n2s}{Sample sizes in control group}
  \item{nsims}{Number of simulations}
  \item{params}{Estimated parameters for simulation from real datasets. Check how to set params in \code{setParams} and \code{estParams}}
  \item{distrib}{Distribution assumption of data generation}
  \item{DAmethod}{Differential abundance analysis methods}
  \item{prevalence}{Non-zero precentage of taxon}
  \item{verbose}{Print out implemented steps if TRUE}
}

\value{
  An object \code{"powmic"} containing:
  \item{pvalues}{A three-dimensional matrix of pvalue (#Taxon by #sample by #simulation)}
  \item{fdrs}{A three-dimensional matrix of fdr (#Taxon by #sample by #simulation)}
  \item{x.abunds}{A three-dimensional matrix of abundance (#Taxon by #sample by #simulation)}
  \item{x.prevs}{A three-dimensional matrix of prevalence (#Taxon by #sample by #simulation)}
  \item{DAids}{A list of vectors. Length of list is #simulation. Each vector contains the number of true Differential Abundant Taxa}
}


\seealso{setParams,assess,summaryAssess,plotStrata,plotStrataAll}
\examples{
\dontrun{
library(powmic)
data(params)
distrib='NB'
setParams=switch(distrib,ZIP=setParams.ZIP, NB=setParams.NB, ZINB=setParams.ZINB)
lmu0=log(params$mu)
lphi0=log(params$phi)
Sigma=params$Sigma
params=setParams(nTaxa=1000,p.DA=0.05,Sigma=Sigma,lmu0=lmu0,lphi0=lphi0,lfc.mu=2,
                        depth=1,sim.seed=123456)
powmic.out=powmic(n1s=c(20,40,60), n2s=c(20,40,60), 
                 params=params,distrib=distrib,DAmethod='edgeR')

assess.out = assess(powmic.out, alpha.type='fdr',alpha.level=0.1,stratify.type='prevalence')
sum.out=summaryAssess(assess.out,assess.type='overall')
}
}





