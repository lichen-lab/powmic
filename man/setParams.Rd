\name{setParams}
\alias{setParams.NB}
\alias{setParams.ZINB}
\alias{setParams.ZIP}
\alias{setParams.DM}
\title{Set parameters for different distributions}
\description{
Four functions for setting parameters for different distributions including Negative Binomial Distribution (NB), Zero-Inflated Negative Binomial Distribution (ZINB), Zero-Inflated Poisson Distribution (ZIP) and Dirichlet-Multinomial Distribution (DM)
}

\usage{
setParams.NB(nTaxa=1000,p.DA=0.05,Sigma=Sigma,lmu0=lmu0,lphi0=lphi0,lfc.mu=2,depth=1,sim.seed=123456)
                        
setParams.ZINB(nTaxa=1000,p.DA=0.05,Sigma=Sigma,lmu0=lmu0,lphi0=lphi0,lp0=lp0,lfc.mu=2,depth=1,sim.seed=123456)

setParams.ZIP(nTaxa=1000,p.DA=0.05,Sigma=Sigma,lmu0=lmu0,lp0=lp0,lfc.mu=2,depth=1,sim.seed=123456) 

setParams.DM(nTaxa=1000,p.DA=0.05,lmu0=lmu0,lphi0=lphi0,lfc.mu=2,libsize=5000,depth=1,sim.seed=123456)
                        
}

\arguments{
  \item{nTaxa}{Number of Taxa in the simulation. If \code{Sigma} is not NULL, nTaxa equals 
  to the dimension of \code{Sigma}}
  \item{p.DA}{Proportion of Differential Abundant (DA) Taxa}
  \item{Sigma}{Estimated microbial correlation matrix when the distribution is "NB", ZINB" or "ZIP"}
  \item{lmu0}{logarithmic mean abundance parameter}
  \item{lphi0}{logarithmic dispersion parameter}
  \item{lp0}{logarithmic zero-inflated probability when the distribution is "ZINB" and "ZIP"}
  \item{lfc.mu}{logarithmic fold change of DA Taxa. Default is 2}
  \item{depth}{Depth factor. depth is multiplied to \code{lmu0} to increase the library size}
  \item{libsize}{Sequence depth when the distribution is "DM"}
  \item{sim.seed}{Simulation seed. Default is 123456}
}

\value{
  A list containing all parameters in the arguments
}


\examples{
\dontrun{
 library(powmic)
 distrib='NB'
 setParams=switch(distrib,ZIP=setParams.ZIP, NB=setParams.NB, ZINB=setParams.ZINB)
 data(params)
 lmu0=log(params$mu)
 lphi0=log(params$phi)
 Sigma=params$Sigma
 sim.params=setParams(Sigma=Sigma,lmu0=lmu0,lphi0=lphi0)
 }
}





