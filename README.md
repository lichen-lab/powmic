# powmic
power assessment in microbiome case-control studies

## Description
Power analysis is essential to decide the sample size of metagenomic sequencing experiments in a case-control study for identifying differentially abundant microbes. However, the complexity of microbiome data characteristics such as excessive zeros, over-dispersion, compositional effect, intrinsically microbial correlations and variable sequencing depths makes the power analysis particularly challenging as the analytical form is usually unavailable. Here, we develop a simulation-based strategy and R package powmic to estimate the empirical statistical power while considering the complexity of data characteristics.

# Maintainer
Li Chen <chen61@iu.edu>


<<<<<<< HEAD
# Install other dependent packages not available on CRAN
```r
#Install R packages for infeerring microbial correlation network
=======
# Install R packages for infeerring microbial correlation network
```r
>>>>>>> Initial commit
#MAGMA
library(devtools)
install_gitlab("arcgl/rmagma")
library(rMAGMA)

#CCREPE
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ccrepe")

#SpiecEasi
library(devtools)
install_github("zdk123/SpiecEasi")
```

<<<<<<< HEAD
=======
# Install other dependent R packages
```{r block2, echo=TRUE,eval=FALSE}
BiocManager::install(c("biomformat","edgeR","DESeq2"))
install.packages(c('ggplot2','gridExtra','lattice','reshape2','MASS','dirmult','nonnest2'))
```


>>>>>>> Initial commit
# Install powmic
```r
install.packages("devtools")
library(devtools)
install_github("lichen-lab/powmic")
```


<<<<<<< HEAD
# Examples

## Estimation module
```r
library(powmic)
#load example dataset
data(params)
params2=estParams(x,Sigma=NULL,method='MAGMA',distrib='NB')
all(params2$Sigma==params$Sigma)
all(params2$mu==params$mu)
all(params2$phi==params$phi)
```

## Simulation module
```r
library(powmic)
#load example dataset
data(params)
#obtain estimated parameters
distrib='NB'
setParams=switch(distrib,ZIP=setParams.ZIP, NB=setParams.NB, ZINB=setParams.ZINB,DM=setParams.DM)
lmu0=log(params$mu)
lphi0=log(params$phi)
Sigma=params$Sigma

#set other parameters
params=setParams(nTaxa=1000,p.DA=0.05,Sigma=Sigma,lmu0=lmu0,lphi0=lphi0,lfc.mu=2,
                 depth=1,sim.seed=123456)
#run simulations
powmic.out=powmic(n1s=c(20,40,60), n2s=c(20,40,60),
                  params=params,distrib=distrib,DAmethod='edgeR')
```

## Assessment module
```r
library(gridExtra)
assess.out = assess(powmic.out, alpha.type='fdr',alpha.level=0.1,stratify.type='prevalence')
sum.out=summaryAssess(assess.out,assess.type='overall')
p1 = plotStrata(assess.out, "power", stratify.by = 'prevalence', 
                is.legend = T)
p2 = plotStrata(assess.out, "fdr", stratify.by = 'prevalence', 
                is.legend = F)
grid.arrange(p1, p2, nrow = 1)
```
=======

# Tutorial

```r
library(powmic)
browseVignettes('powmic')
```

**Tutorial** is available [**here**](http://htmlpreview.github.io/?https://github.com/lichen-lab/powmic/blob/master/inst/doc/powmic.html)


>>>>>>> Initial commit






