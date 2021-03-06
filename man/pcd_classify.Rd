% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcd.R
\name{pcd_classify}
\alias{pcd_classify}
\title{Classify an unlabelled data based on a PCD of standard cover. The class cover is
composed of only simplical proximity regions.}
\usage{
pcd_classify(data, graph_pcd)
}
\arguments{
\item{data}{An m-by-d matrix of the test data set.}

\item{graph_pcd}{A proximity catch digraph (PCD).}
}
\value{
Predicted labels of the unlabelled data set.
}
\description{
Classify an unlabelled data based on a PCD of standard cover. The class cover is
composed of only simplical proximity regions.
}
\examples{

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
}
