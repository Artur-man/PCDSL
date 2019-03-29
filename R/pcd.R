#' The PCD classifier based on the cover of the two (or more) classes wherein the cover of
#' the target class is restricted to the inside of the convex hull of the non-target class.
#'
#' @param data    An n-by-d matrix of the training data set.
#' @param classes A vector of length n indicating the labels of the classes.
#' @param map     The Proximity Map associated with the classifier. "pe" for the Proportional-edge proximity maps.
#'                and "cs" for the Central-similarity proximity maps.
#' @param param   The value of the parameter associated with the proximity map.
#'
#' @return        A proximity catch digraph (PCD)
#'
#' @examples
#'
#' # input parameters
#' ntest <- 100     # test data size for each class
#' nx <- 300        # training data size of x class (majority)
#' r <- 0.1         # Imbalance Ratio
#' de <- 0.5        # delta, the overlapping parameter
#' dimx <- 2        # number of dimensions
#'
#' # training the classifier
#' set.seed(1)
#' x0 <- matrix(runif(dimx*nx,0,1),nrow=nx)
#' x1 <- matrix(runif(dimx*nx*r,de,1+de),nrow=nx*r)
#' x <- rbind(x0,x1)
#' classes <- rep(1:2,c(nx,nx*r))
#' graph_pcd <- pcd_inhull_classifier(x,classes,map="pe",param = 1)
pcd_inhull_classifier <- function(data, classes, map = "pe", param = 1)
{
  cls <- unique(classes)
  nc <- length(cls)
  Simp <- D <- G <- C <- list()
  length(G) <- length(C) <- nc
  length(Simp) <- length(D) <- nc

  for (i in 1:nc)
  {
    z <- classes == cls[i]
    x <- data[z, ]
    y <- data[!z, ]

    Te <- delaunayn(y)
    Ts <- tsearchn(y, Te, x)
    indx <- which(!is.na(Ts$idx))
    indT <- Ts$idx[indx]
    baryx <- Ts$p[indx, ]
    if (!is.matrix(baryx)) baryx <- matrix(baryx, nrow = 1)

    for (j in unique(indT))
    {
      ind <- which(indT == j)
      xT <- x[indx[ind], ]
      yT <- y[Te[j, ], ]
      baryT <- baryx[ind, ]
      if (map == "pe") pcd_simp <- pcd_pe_simp(xT, yT, baryT, r = param)
      if (map == "cs") pcd_simp <- pcd_cs_simp(xT, yT, baryT, tau = param)

      nD <- indx[ind[pcd_simp$D]]
      nSimp <- pcd_simp$Simp
      D[[i]] <- append(D[[i]], nD)
      Simp[[i]] <- c(Simp[[i]], nSimp)
    }

    C[[i]] <- which(classes == cls[i])
  }

  return(list(D = D, C = C, Simp = Simp, classes = cls))
}

#' Classify an unlabelled data set with a PCD wherein the cover of
#' the target class is restricted to inside of the convex hull of the non-target class.
#' Returns NA for points outside of the convex hull of the non-target classes.
#'
#' @param data An m-by-d matrix of the test data set.
#' @param graph A proximity catch digraph (PCD).
#'
#' @return Predicted labels of the test data set.
#'
#' @examples
#'
#' # input parameters
#' ntest <- 100     # test data size for each class
#' nx <- 300        # training data size of x class (majority)
#' r <- 0.1         # Imbalance Ratio
#' de <- 0.5        # delta, the overlapping parameter
#' dimx <- 2        # number of dimensions
#'
#' # training the classifier
#' set.seed(1)
#' x0 <- matrix(runif(dimx*nx,0,1),nrow=nx)
#' x1 <- matrix(runif(dimx*nx*r,de,1+de),nrow=nx*r)
#' x <- rbind(x0,x1)
#' classes <- rep(1:2,c(nx,nx*r))
#' graph_pcd <- pcd_inhull_classifier(x,classes,map="pe",param = 1)
#'
#' # testing
#' tx0 <- matrix(runif(dimx*ntest,0,1),nrow=ntest)
#' tx1 <- matrix(runif(dimx*ntest,de,1+de),nrow=ntest)
#' tx <- rbind(tx0,tx1)
#' tclsdata <- rep(1:2,rep(ntest,2))
#' predicted_pcd_tx <- pcd_inhull_classify(tx,graph_pcd)
#'
#' # classifying with the alternative classifier for those points that are not
#' # classified
#' library(e1071)
#' svm_x <- svm(x = x, y = factor(classes), kernel="radial", gamma = 0.1)
#' predicted_svm_tx <- predict(svm_x,tx)
#' predicted_pcd_tx[is.na(predicted_pcd_tx)] <-predicted_svm_tx[is.na(predicted_pcd_tx)]
pcd_inhull_classify <- function(data, graph)
{
  ncn <- ncol(data)
  cls <- graph$classes
  dist_matrix <- NULL
  D <- list()
  length(D) <- length(cls)

  for(i in 1:length(cls)) {
    temp <- graph$D[[i]]
    if(!is.null(temp)) D[[i]] <- temp
    else D[[i]] <- numeric(0)
  }
  nc <- unlist(lapply(D,length))
  if(!is.null(nc)) cl <- rep(cls,nc)

  result <- rep(0,nrow(data))

  for(i in 1:length(cls))
  {
    Simp <- graph$Simp[[i]]
    if(length(Simp)>0){
      temp <- sapply(Simp,function(simp){
        baryx <- cart2bary(simp,data)
        smallbary <- apply(baryx,1,function(t){
          if(any(t < 0)){
            return(Inf)
          } else{
            return(1-(ncn+1)*min(t))
          }
        })
      }, simplify=TRUE)
      dist_matrix <- cbind(dist_matrix,temp)
    }
  }

  if(!is.null(dist_matrix)){
    dist_matrix <- matrix(dist_matrix,nrow=nrow(data))
    in_matrix <- apply(dist_matrix,1,which.min)
    dist_res <- apply(dist_matrix,1,min)
    in_dist <- which(dist_res!=Inf)

    if(length(in_dist) > 0){
      if(length(in_dist)==1) dist_matrix <- matrix(dist_matrix[in_dist,],nrow=1)
      else dist_matrix <- matrix(dist_matrix[in_dist,],nrow=length(in_dist))

      result[in_dist] <- apply(dist_matrix,1,function(z){
        return(cl[which.min(z)])
      })

      result[-in_dist] <- NA

    } else {
      result <- rep(NA,length(cls))
    }
  } else {
    result <- rep(NA,length(cls))
  }

  return(result)
}

#' The PCD classifier based on the cover of the two (or more) classes wherein the
#' points inside the convex hull are classified with PCDs, and ones that are outside
#' with CCCDs
#'
#' @param data An m-by-d matrix of the training data set.
#' @param classes A vector of length n indicating the labels of the classes.
#' @param map The Proximity Map associated with the classifier. "pe" for the Proportional-edge proximity maps
#'               and "cs" for the Central-similarity proximity maps.
#' @param p_pcd  The value of the parameter associated with the proximity map of the PCD.
#' @param p_cccd  The value of the parameter associated with the CCCD.
#'
#' @return A proximity catch digraph and a class cover catch digraph with the minimum dominating set and
#'         the associated proximity regions.
#'
#' @examples
#'
#' # input parameters
#' ntest <- 100     # test data size for each class
#' nx <- 300        # training data size of x class (majority)
#' r <- 0.1         # Imbalance Ratio
#' de <- 0.5        # delta, the overlapping parameter
#' dimx <- 2        # number of dimensions
#'
#' # training the classifier
#' set.seed(1)
#' x0 <- matrix(runif(dimx*nx,0,1),nrow=nx)
#' x1 <- matrix(runif(dimx*nx*r,de,1+de),nrow=nx*r)
#' x <- rbind(x0,x1)
#' classes <- rep(1:2,c(nx,nx*r))
#' graph_pcd <- pcd_cccd_classifier(x,classes,map="pe",p_pcd=1,p_cccd=1)
pcd_cccd_classifier <- function(data,classes,map="pe",p_pcd=1,p_cccd=1)
{
  cls <- unique(classes)
  nc <- length(cls)
  D_cccd <- G <- C <- Ti <- R <- Simp <- D_pcd <- list()
  length(G) <- length(C) <- length(Ti) <- length(Simp) <- nc
  length(D_cccd) <- length(R) <- length(D_pcd) <- nc

  for (i in 1:nc)
  {
    z <- classes == cls[i]
    C[[i]] <- which(classes==cls[i])
    x <- data[z,]
    y <- data[!z,]

    Te <- delaunayn(y)
    Ti[[i]] <- Te
    Ts <- tsearchn(y, Te, x)
    indx <- which(!is.na(Ts$idx))
    indnx <- which(is.na(Ts$idx))
    indT <- Ts$idx[indx]
    baryx <- Ts$p[indx,]
    if(!is.matrix(baryx)) baryx <- matrix(baryx,nrow=1)

    for(j in unique(indT))
    {
      ind <- which(indT==j)
      xT <- x[indx[ind],]
      yT <- y[Te[j,],]
      baryT <- baryx[ind,]
      if(map=="pe") pcd_simp <- pcd_pe_simp(xT,yT,baryT,r=p_pcd)
      if(map=="cs") pcd_simp <- pcd_cs_simp(xT,yT,baryT,tau=p_pcd)

      nD    <- indx[ind[pcd_simp$D]]
      D_pcd[[i]] <- append(D_pcd[[i]],nD)
      nSimp  <- pcd_simp$Simp
      Simp[[i]] <- c(Simp[[i]],nSimp)
    }
    if(length(indnx) > 0)
    {
      G_cccd <- cccd(x[indnx,], y=y)
      Dom <- dominate(G_cccd)
      D_cccd[[i]] <- C[[i]][indnx[Dom]]
      R[[i]] <- G_cccd$R[Dom]
    }
  }

  return(list(D_pcd=D_pcd,D_cccd=D_cccd,C=C,p_pcd=p_pcd,p_cccd=p_cccd,
              Simp=Simp,R=R,Ti=Ti,classes=cls))
}

#' Classify an unlabelled data set with a PCD of composite cover. The class cover is
#' composed of simplical and spherical proximity regions.
#'
#' @param data An m-by-d matrix of the test data set.
#' @param dvalid The distance matrix for the distances between the test and training data.
#' @param graph_pcd A proximity catch digraph (PCD).
#'
#' @return Predicted labels of the test data set.
#'
#' @examples
#'
#' # input parameters
#' ntest <- 100     # test data size for each class
#' nx <- 300        # training data size of x class (majority)
#' r <- 0.1         # Imbalance Ratio
#' de <- 0.5        # delta, the overlapping parameter
#' dimx <- 2        # number of dimensions
#'
#' # training the classifier
#' set.seed(1)
#' x0 <- matrix(runif(dimx*nx,0,1),nrow=nx)
#' x1 <- matrix(runif(dimx*nx*r,de,1+de),nrow=nx*r)
#' x <- rbind(x0,x1)
#' classes <- rep(1:2,c(nx,nx*r))
#' graph_pcd <- pcd_cccd_classifier(x,classes,map="pe",p_pcd=1,p_cccd=1)
#'
#' # testing
#' tx0 <- matrix(runif(dimx*ntest,0,1),nrow=ntest)
#' tx1 <- matrix(runif(dimx*ntest,de,1+de),nrow=ntest)
#' tx <- rbind(tx0,tx1)
#' tclsdata <- rep(1:2,rep(ntest,2))
#' library(flexclust)
#' d_tx_x <- as.matrix(dist2(tx,x))
#' predicted_pcd_tx <- pcd_cccd_classify(tx,d_tx_x,graph_pcd)
pcd_cccd_classify <- function(data,dvalid,graph_pcd)
{
  ncn <- ncol(data)
  cls <- graph_pcd$classes
  dist_matrix <- NULL

  D_pcd <- list()
  for(i in 1:length(cls)) {
    temp <- graph_pcd$D_pcd[[i]]
    if(!is.null(temp)) D_pcd[[i]] <- temp
    else D_pcd[[i]] <- numeric(0)
  }
  ncp <- unlist(lapply(D_pcd,length))
  if(!is.null(ncp)) cl_pcd <- rep(cls,ncp)

  D_cccd <- graph_pcd$D_cccd
  nc <- unlist(lapply(D_cccd,length))
  cl_cccd <- rep(cls,nc)
  R <- unlist(graph_pcd$R)
  D_cccd <- unlist(D_cccd)
  distd <- dvalid[,D_cccd]

  if(!is.null(ncp)) cl <- c(cl_pcd,cl_cccd)
  else cl <- cl_cccd

  for(i in 1:length(cls))
  {
    Simp <- graph_pcd$Simp[[i]]
    if(length(Simp)>0){
      temp <- sapply(Simp,function(simp){
        baryx <- cart2bary(simp,data)
        smallbary <- apply(baryx,1,function(t){
          return(1-(ncn+1)*min(t))
        })
      }, simplify=TRUE)
      dist_matrix <- cbind(dist_matrix,temp)
    }
  }

  dist.cccd <- t(apply(distd,1,function(z){
    return(z/R)
  }))
  dist_matrix <- cbind(dist_matrix,dist.cccd)
  result <- apply(dist_matrix,1,function(z){
    return(cl[which.min(z)])
  })

  return(result)
}

#' The PCD classifier based on the cover of the two (or more) classes wherein the target class points inside the
#' convex hull are classified with n-simplices, and the ones that are outside with outer simplices.
#'
#' @param data An n-by-d matrix of the training data set.
#' @param classes A vector of length n indicating the labels of the classes.
#' @param map The Proximity Map associated with the classifier. "pe" for the Proportional-edge proximity maps
#'               and "cs" for the Central-similarity proximity maps.
#' @param p_pcd  The value of the parameter associated with the proximity map.
#'
#' @return A proximity catch digraph (PCD).
#'
#' @examples
#'
#' # input parameters
#' ntest <- 100     # test data size for each class
#' nx <- 300        # training data size of x class (majority)
#' r <- 0.1         # Imbalance Ratio
#' de <- 0.5        # delta, the overlapping parameter
#' dimx <- 2        # number of dimensions
#'
#' # training the classifier
#' set.seed(1)
#' x0 <- matrix(runif(dimx*nx,0,1),nrow=nx)
#' x1 <- matrix(runif(dimx*nx*r,de,1+de),nrow=nx*r)
#' x <- rbind(x0,x1)
#' classes <- rep(1:2,c(nx,nx*r))
#' graph_pcd <- pcd_classifier(x,classes,map="pe",p_pcd=1)
pcd_classifier <- function(data,classes,map="pe",p_pcd=1)
{
  cls <- unique(classes)
  nc <- length(cls)
  G <- C <- Ti <- D.in <- D.out <- Simp <- Osimp <- list()
  length(G) <- length(C) <- length(Ti) <- nc
  length(Simp) <- length(Osimp) <- nc
  length(D.in) <- length(D.out) <- nc

  for (i in 1:nc)
  {
    z <- classes == cls[i]
    C[[i]] <- which(classes==cls[i])
    x <- data[z,]
    y <- data[!z,]

    Te <- delaunayn(y)
    Ti[[i]] <- Te
    Ts <- tsearchn(y, Te, x)
    indx <- which(!is.na(Ts$idx))
    indnx <- which(is.na(Ts$idx))
    indT <- Ts$idx[indx]
    baryx <- Ts$p[indx,]
    if(!is.matrix(baryx)) baryx <- matrix(baryx,nrow=1)

    for(j in unique(indT))
    {
      ind <- which(indT==j)
      xT <- x[indx[ind],]
      yT <- y[Te[j,],]
      baryT <- baryx[ind,]

      if(map=="pe") pcd_simp <- pcd_pe_simp(xT,yT,baryT,r=p_pcd)
      if(map=="cs") pcd_simp <- pcd_cs_simp(xT,yT,baryT,tau=p_pcd)

      nD    <- indx[ind[pcd_simp$D]]
      D.in[[i]] <- append(D.in[[i]],nD)
      nSimp  <- pcd_simp$Simp
      Simp[[i]] <- c(Simp[[i]],nSimp)
    }

    if(length(indnx) > 0)
    {
      ETe <- outer_delaunayn(y)
      indTe <- etsearchn(y,ETe$Te,x[indnx,])

      for(k in 1:length(p_pcd))
      {
        for(j in unique(indTe))
        {
          ind <- which(indTe==j)
          xT <- x[indnx[ind],]
          yT <- ETe$Te[[j]]
          CM <- ETe$CM

          if(map=="pe") pcd_Osimp <- pcd_pe_Osimp(xT,yT,r=p_pcd)
          if(map=="cs") pcd_Osimp <- pcd_cs_Osimp(xT,yT,CM,tau=p_pcd)

          nD    <- indnx[ind[pcd_Osimp$D]]
          D.out[[i]] <- append(D.out[[i]],nD)
          nSimp  <- pcd_Osimp$Simp
          Osimp[[i]] <- append(Osimp[[i]],nSimp)
        }
      }
    }
  }

  return(list(D.in=D.in,D.out=D.out,C=C,p_pcd=p_pcd,Simp=Simp,Osimp=Osimp,Ti=Ti,classes=cls))
}

#' Classify an unlabelled data based on a PCD of standard cover. The class cover is
#' composed of only simplical proximity regions.
#'
#' @param data An m-by-d matrix of the test data set.
#' @param graph_pcd A proximity catch digraph (PCD).
#'
#' @return Predicted labels of the unlabelled data set.
#'
#' @examples
#'
#' # input parameters
#' ntest <- 100     # test data size for each class
#' nx <- 300        # training data size of x class (majority)
#' r <- 0.1         # Imbalance Ratio
#' de <- 0.5        # delta, the overlapping parameter
#' dimx <- 2        # number of dimensions
#'
#' # training the classifier
#' set.seed(1)
#' x0 <- matrix(runif(dimx*nx,0,1),nrow=nx)
#' x1 <- matrix(runif(dimx*nx*r,de,1+de),nrow=nx*r)
#' x <- rbind(x0,x1)
#' classes <- rep(1:2,c(nx,nx*r))
#' graph_pcd <- pcd_classifier(x,classes,map="pe",p_pcd=1)
#'
#' # testing
#' tx0 <- matrix(runif(dimx*ntest,0,1),nrow=ntest)
#' tx1 <- matrix(runif(dimx*ntest,de,1+de),nrow=ntest)
#' tx <- rbind(tx0,tx1)
#' tclsdata <- rep(1:2,rep(ntest,2))
#' predicted_pcd_tx <- pcd_classify(tx,graph_pcd)
pcd_classify <- function(data,graph_pcd)
{
  ncn <- ncol(data)
  cls <- graph_pcd$classes
  dist_matrix <- NULL

  D <- list()
  for(i in 1:length(cls))
    D[[i]] <- c(graph_pcd$D.in[[i]],graph_pcd$D.out[[i]])
  nc <- unlist(lapply(D,length))
  cl <- rep(cls,nc)

  for(i in 1:length(cls))
  {
    Simp <- graph_pcd$Simp[[i]]
    if(length(Simp)>0){
      temp <- sapply(Simp,function(simp){
        baryx <- cart2bary(simp,data)
        if(!is.matrix(baryx)) baryx <- matrix(baryx,nrow=1)
        smallbary <- apply(baryx,1,function(t){
          return(1-(ncn+1)*min(t))
        })
      }, simplify=TRUE)
      dist_matrix <- cbind(dist_matrix,temp)
    }

    Osimp <- graph_pcd$Osimp[[i]]
    if(length(Osimp) > 0){
      temp <- sapply(Osimp,getDisttoOsimp,x=data)
      dist_matrix <- cbind(dist_matrix,temp)
    }
  }

  result <- apply(dist_matrix,1,function(z){
    return(cl[which.min(z)])
  })

  return(result)
}
