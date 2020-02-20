\name{assess}
\alias{assess}
\title{Perform power assessment}
\description{
    Assess the simulation results from powmic by providing overall and stratified FDR (False Discovery Rate), TP (True Positive), FP (False Positive),  TPR (True Positive Rate), FP (False Positive Rate)
}

\usage{
assess(powmic.out, alpha.type=c("fdr","pval"), alpha.level=0.1,stratify.type=c("prevalence","abundance"))
}

\arguments{
  \item{powmic.out}{An \code{powmic} object}
  \item{alpha.type}{Power assessment using FDR ("fdr") or pvalue ("pval")}
  \item{alpha.level}{Significant level. Default is 0.1}
  \item{stratify.type}{OTU is stratified by either "prevalence" or "abundance"}
}



\seealso{summaryAssess, plotStrata, plotStrataAll}
\examples{
 powmic.out=powmic(n1s=c(20,40,60), n2s=c(20,40,60), nsims=10, 
                 sim.params=sim.params,distrib=distrib,DAmethod='edgeR')
 assess.out = assess(powmic.out, alpha.type='fdr',alpha.level=0.1,stratify.type='prevalence')
}






