\name{summaryAssess}
\alias{summaryAssess}
\title{Summarize the power assessment criteria}
\description{
    Summarize power assessment criteria including overall or stratified FDR (False Discovery Rate), TP (True Positive), FP (False Positive),  TPR (True Positive Rate), FP (False Positive Rate)
}

\usage{
summaryAssess(assess.out,assess.type=c('overall','stratified'))
}

\arguments{
  \item{assess.out}{Output object from assess}
  \item{assess.type}{'overall' or 'stratified'}
}



\seealso{assess}
\examples{
 \dontrun{
assess.out = assess(powmeta.out, alpha.type='fdr',alpha.level=0.1,stratify.type='prevalence')
sum.out=summaryAssess(assess.out,assess.type='overall')
}

}






