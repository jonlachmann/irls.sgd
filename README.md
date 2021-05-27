# irls.sgd
Approximately fit a GLM model using Subsampling IRLS and SGD

The `irls.sgd` package provides functions to acquire an approximate estimate of the MLE, specifically tailored to be used in Laplace approximations for BGNLM in scenarios when large amounts of data is available.

# Installation and getting started
To install and load the package, just run
```
library(devtools)
install_github("jonlachmann/irls.sgd", force=T, build_vignettes=T)
library(irls.sgd)
```
With the package loaded, functions of interest are `irls.sgd`, `sgd` and `glm.sgd`. Documentation is available once installed by using the built in help in R.

