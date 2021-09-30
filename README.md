# evolCCM
R package for CCM (Community Coevolution Model).

## Introduction
This is the R implementation for CCM (Community Coevolution Model) from the manuscript "The Community Coevolution Model with Application to the Study of Evolutionary Relationships between Genes based on Phylogenetic Profiles". This package was written in R v4.0.2 and includes the main functions for profiles simulation, profiles visualization and CCM estimation.

## Dependencies

`evolCCM` requires the following R packages:

- [ape](https://cran.r-project.org/web/packages/ape/)

- [gplots](https://cran.r-project.org/web/packages/gplots/index.html)

```
# install.packages("ape")
# install.packages("gplots")
library(ape)
library(gplots)
```


## Installation

Install the package `evolCCM` from github using `devtools`:

```
# install.packages("devtools")
library(devtools)
devtools::install_github("beiko-lab/evolCCM")
library(evolCCM)
```


## Documentation

For detailed instructions, please see the package [manual](https://github.com/beiko-lab/evolCCM/blob/main/evolCCM_manual.pdf). 

