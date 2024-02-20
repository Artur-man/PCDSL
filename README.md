
<!-- README.md is generated from README.Rmd. Please edit that file -->
PCDSL
=====

PCDSL: The Proximity Catch Digraphs for Statistical Learning

Artür Manukyan, Elvan Ceyhan

Classification methods based on Proximity Catch Digraphs. Proximity Catch Digraphs (PCDs) are special types of proximity graphs.

Webpage
------------

http://webhome.auburn.edu/~ezc0066/PCDwebpage/index.html


Installation
------------

``` r
library(devtools)
devtools::install_github("Artur-man/PCDSL")
```

Example
-------

``` r
# input parameters
ntest <- 100     # test data size for each class
nx <- 300        # training data size of x class (majority)
r <- 0.1         # Imbalance Ratio
de <- 0.5        # delta, the overlapping parameter
dimx <- 2        # number of dimensions

# training the classifier
set.seed(1)
x0 <- matrix(runif(dimx*nx,0,1),nrow=nx)
x1 <- matrix(runif(dimx*nx*r,de,1+de),nrow=nx*r)
x <- rbind(x0,x1)
classes <- rep(1:2,c(nx,nx*r))
graph_pcd <- pcd_classifier(x,classes,map="pe",p_pcd=1)

# testing
tx0 <- matrix(runif(dimx*ntest,0,1),nrow=ntest)
tx1 <- matrix(runif(dimx*ntest,de,1+de),nrow=ntest)
tx <- rbind(tx0,tx1)
tclsdata <- rep(1:2,rep(ntest,2))
predicted_pcd_tx <- pcd_classify(tx,graph_pcd)
```
