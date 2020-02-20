\name{estParams}
\alias{estParams}
\title{estimate paramters from real datasets based on different distribution assumptions}
\description{
    Estimate paramters from real datasets based on different distribution assumptions including NB (Negative Binomial), ZINB (Zero-Inflated Negative Binomial), ZIP (Zero-Inflated Poisson) and Dirichlet-Multinomial Distribution (DM). The microbial correlation network can be provided or estimated from different methods including 'MAGMA','CCLasso','SparCC','SPIEC-EASI' and 'CCREPE'.
}

\usage{
estParams(otu_table,Sigma=NULL,
	method=c('MAGMA','CCLasso','SparCC','SPIECEASI','CCREPE'),
	distrib=c('NB','ZINB','ZIP','DM'))
}

\arguments{
  \item{x}{taxon abundance matrix with rows as samples and columns as taxa}
  \item{Sigma}{correlation matrix among taxa}
  \item{method}{methods for inferring Sigma if not provided}
  \item{distrib}{distribution can be specified as 'NB','ZINB','ZIP' or 'DM'}
}

\value{
  A list containing all parameters: dispersion (\code{phi}) and mean (\code{mu}) for NB; dispersion (\code{phi}), mean (\code{mu}) and zero inflation probability (\code{p0}) for ZINB; mean (\code{mu}) and zero inflation probability (\code{p0}) for ZIP; mean (\code{mu}) and dispersion (\code{phi}) for DM
}


\examples{
 \dontrun{
 library(powmic)
 x=params$otu.mat
 params=estParams(x,Sigma=NULL,method='CCLasso')
 library(rMAGMA)
 data(params)
 params=estParams(x,Sigma=NULL,method='MAGMA',distrib='NB')
 library(SpiecEasi)
 params=estParams(x,Sigma=NULL,method='SPIECEASI')
 params=estParams(x,Sigma=NULL,method='SparCC')
 library(ccrepe)
 params=estParams(x,Sigma=NULL,method='CCREPE')
 }
}











