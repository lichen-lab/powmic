# powmic
power assessment in microbiome case-control studies

## Description
Power analysis is essential to decide the sample size of metagenomic sequencing experiments in a case-control study for identifying differentially abundant microbes. However, the complexity of microbiome data characteristics such as excessive zeros, over-dispersion, compositional effect, intrinsically microbial correlations and variable sequencing depths makes the power analysis particularly challenging as the analytical form is usually unavailable. Here, we develop a simulation-based strategy and R package powmic to estimate the empirical statistical power while considering the complexity of data characteristics.

# Maintainer
Li Chen <chen61@iu.edu>


# Install R packages for infeerring microbial correlation network
```r
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

# Install other dependent R packages
```{r block2, echo=TRUE,eval=FALSE}
BiocManager::install(c("biomformat","edgeR","DESeq2"))
install.packages(c('ggplot2','gridExtra','lattice','reshape2','MASS','dirmult','nonnest2'))
```

# Install powmic
```r
install.packages("devtools")
library(devtools)
install_github("lichen-lab/powmic")
```


# Tutorial

```r
library(powmic)
browseVignettes('powmic')
```

**Tutorial** is available [**here**](http://htmlpreview.github.io/?https://github.com/lichen-lab/powmic/blob/master/inst/doc/powmic.html)








