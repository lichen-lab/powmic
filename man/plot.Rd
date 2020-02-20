\name{plot}
\alias{plotStrata}
\alias{plotStrataAll}
\title{Plot the stratified power assessment criteria}
\description{
    Plot the stratified Power, FDR (False Discovery Rate), TP (True Positive), FP (False Positive) based on prevalence or abundance.
}

\usage{
plotStrata(assess.out, figure.type=c('power','fdr','type1error','tp','fp','otu'),stratify.by=c('prevalence','abundance'),is.errorbar=TRUE,is.legend=TRUE) 
plotStrataAll(assess.out, figure.type=c('power','fdr','type1error','tp','fp','otu'),stratify.by=c('prevalence','abundance'),is.errorbar=TRUE,is.legend=TRUE)
}

\arguments{
  \item{assess.out}{An object output by \code{assess} function}
  \item{figure.type}{Figure type could be Power("power"), FDR ("fdr"), TypeI error ("type1error"), TP("tp"), FP("tp"), #OTU ("otu")}
  \item{stratify.type}{Taxon is stratified by either "prevalence" or "abundance"}
  \item{is.errorbar}{Standard error of all simulations is plot}
  \item{is.legend}{Legend of sample sizes is plot}
}


\examples{
\dontrun{
assess.out = assess(powmeta.out, alpha.type="fdr", alpha.level=0.1,stratify.type='prevalence')
plotStrataAll(assess.out,stratify.by='prevalence')
plotStrata(assess.out,figure.type='power',stratify.by='prevalence')
 }
}

\keyword{plot}







