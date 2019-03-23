#' The PCD classifier based on the proximity regions of the target class which are restricted to the inside of the convex hull of
#' the non-target class
#'
#' @param datax The labelled data set
#' @param classes A categorical variable indicating the classes of the points in the data set
#' @param method The Proximity Map associated with the classifier. "pe" for the Proportional-edge proximity maps
#'               and "cs" for the Central-similarity proximity maps.
#' @param param  The value of the parameter associated with the proximity map.
#'
#' @return A proximity catch digraph with the minimum dominating set and the associated proximity regions
pcd_inhull_classifier <- function(datax, classes, method = "pe", param = 1)
{
  cls <- unique(classes)
  nc <- length(cls)
  Simp <- D <- G <- C <- list()
  length(G) <- length(C) <- nc
  length(Simp) <- length(C) <- nc

  for (i in 1:nc)
  {
    z <- classes == cls[i]
    x <- datax[z, ]
    y <- datax[!z, ]

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
      if (method == "pe") pcd_simp <- pcd_pe_simp(xT, yT, baryT, r = param)
      if (method == "cs") pcd_simp <- pcd_cs_simp(xT, yT, baryT, tau = param)

      nD <- indx[ind[pcd_simp$D]]
      nSimp <- pcd_simp$Simp
      D[[i]] <- append(D[[i]], nD)
      Simp[[i]] <- c(Simp[[i]], nSimp)
    }

    C[[i]] <- which(classes == cls[i])
  }

  return(list(D = D, C = C, P = param, Simp = Simp, classes = cls))
}

#' Classify a set of unlabelled data based on a PCD. Returns NA for points outside of the convex hull of the non-target classes.
#'
#' @param datan The unlabelled data set
#' @param graph A proximity catch digraph (PCD)
#' @param param The value of the parameter of the associated proximity map of the PCD
#'
#' @return Predicted labels of the unlabelled data set.
pcd_inhull_classify <- function(datan, graph, param)
{
  ncn <- ncol(datan)
  cls <- graph$classes
  dist_matrix <- NULL

  par <- which(graph$P==param)
  for(i in 1:length(cls)) {
    temp <- graph$D[[i]]
    if(!is.null(temp)) D[[i]] <- temp
    else D[[i]] <- numeric(0)
  }
  nc <- unlist(lapply(D,length))
  if(!is.null(nc)) cl <- rep(cls,nc)

  result <- rep(0,nrow(datan))

  for(i in 1:length(cls))
  {
    Simp <- graph$Simp[[i,par]]
    if(length(Simp)>0){
      temp <- sapply(Simp,function(simp){
        baryx <- cart2bary(simp,datan)
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
    dist_matrix <- matrix(dist_matrix,nrow=nrow(datan))
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

#' The PCD classifier based on the proximity regions of the target class where points inside the
#' convex hull are classified with PCDs, and outside with CCCDs
#'
#' @param datax The labelled data set
#' @param classes A categorical variable indicating the classes of the points in the data set
#' @param method The Proximity Map associated with the classifier. "pe" for the Proportional-edge proximity maps
#'               and "cs" for the Central-similarity proximity maps.
#' @param p_pcd  A set of values of the parameter associated with the proximity map of the PCD.
#' @param p_cccd  The value of the parameter associated with the CCCD.
#'
#' @return A proximity catch digraph and a class cover catch digraph with the minimum dominating set and
#'         the associated proximity regions.
pcd_cccd_multiclass_classifier <- function(datax,classes,method="pe",p_pcd=1,p_cccd=1)
{
  ddatax <- as.matrix(dist(datax))
  cls <- unique(classes)
  nc <- length(cls)
  D_cccd <- G <- C <- Ti <- R <- list()
  length(G) <- length(C) <- length(Ti) <- length(D_cccd) <- length(R) <- nc
  Simp <- matrix(list(), nrow=nc, ncol=length(p_pcd))
  D_pcd <- matrix(list(), nrow=nc, ncol=length(p_pcd))

  for (i in 1:nc)
  {
    z <- classes == cls[i]
    C[[i]] <- which(classes==cls[i])
    x <- datax[z,]
    y <- datax[!z,]
    dx <- ddatax[z,z]
    dy <- ddatax[!z,z]

    # Delaunayn triangulation info and points
    Te <- delaunayn(y)
    Ti[[i]] <- Te
    Ts <- tsearchn(y, Te, x)
    indx <- which(!is.na(Ts$idx))
    indnx <- which(is.na(Ts$idx))
    indT <- Ts$idx[indx]
    baryx <- Ts$p[indx,]
    if(!is.matrix(baryx)) baryx <- matrix(baryx,nrow=1) # if a single point, create matrix

    # PCDs
    for(k in 1:length(p_pcd))
    {
      # for all triangles
      for(j in unique(indT))
      {
        # which x points are dominating points in the triangles
        ind <- which(indT==j)
        xT <- x[indx[ind],]
        yT <- y[Te[j,],]
        baryT <- baryx[ind,]
        if(method=="pe") pcd_simp <- pcd_pe_simp(xT,yT,baryT,r=p_pcd[k])
        if(method=="cs") pcd_simp <- pcd_cs_simp(xT,yT,baryT,tau=p_pcd[k])

        # add new dominating points, their triangle info, their vertex region
        nD    <- indx[ind[pcd_simp$D]]
        D_pcd[[i,k]] <- append(D_pcd[[i,k]],nD)
        nSimp  <- pcd_simp$Tri
        Simp[[i,k]] <- c(Simp[[i,k]],nSimp)
      }
    }

    # CCCDs
    if(length(indnx) > 0)
    {
      G_cccd <- cccd(dxx=dx[indnx,indnx], dyx=dy[,indnx])
      Dom <- dominate(G_cccd)
      D_cccd[[i]] <- C[[i]][indnx[Dom]]
      R[[i]] <- G_cccd$R[Dom]
    }
  }

  return(list(D_pcd=D_pcd,D_cccd=D_cccd,C=C,p_pcd=p_pcd,p_cccd=p_cccd,Simp=Simp,R=R,Ti=Ti,classes=cls))
}

#' Classify a set of unlabelled data based on a PCD and a CCCD.
#'
#' @param datan The unlabelled data set
#' @param dvalid The distance matrix from the dominating set of CCCD to points in datan
#' @param graph_pcd A proximity catch digraph (PCD)
#' @param p_pcd The value of the parameter of the associated proximity map of the PCD
#'
#' @return Predicted labels of the unlabelled data set.
pcd_cccd_multiclass_classify <- function(datan,dvalid,graph_pcd,p_pcd)
{
  ncn <- ncol(datan)
  cls <- graph_pcd$classes
  dist_matrix <- NULL
  par <- which(graph_pcd$p_pcd==p_pcd)

  D_pcd <- list()
  for(i in 1:length(cls)) {
    temp <- graph_pcd$D.pcd[[i,par]]
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
    Simp <- graph_pcd$Simp[[i,par]]
    if(length(Simp)>0){
      temp <- sapply(Simp,function(simp){
        baryx <- cart2bary(simp,datan)
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

#' The PCD classifier based on the proximity regions of the target class where points inside the
#' convex hull are classified with n-simplicies, and outside with outer simplicies
#'
#' @param datax The labelled data set
#' @param classes A categorical variable indicating the classes of the points in the data set
#' @param method The Proximity Map associated with the classifier. "pe" for the Proportional-edge proximity maps
#'               and "cs" for the Central-similarity proximity maps.
#' @param p_pcd  A set of values of the parameter associated with the proximity map of the PCD.
#'
#' @return A proximity catch digraph with the minimum dominating set and
#'         the associated proximity regions.
pcd_multiclass_classifier <- function(datax,classes,method="pe",p_pcd=1)
{
  cls <- unique(classes)
  nc <- length(cls)
  G <- C <- Ti <- list()
  length(G) <- length(C) <- length(Ti) <- nc
  Simp <- matrix(list(), nrow=nc, ncol=length(p_pcd))
  Osimp <- matrix(list(), nrow=nc, ncol=length(p_pcd))
  D.in <- matrix(list(), nrow=nc, ncol=length(p_pcd))
  D.out <- matrix(list(), nrow=nc, ncol=length(p_pcd))

  for (i in 1:nc)
  {
    z <- classes == cls[i]
    C[[i]] <- which(classes==cls[i])
    x <- datax[z,]
    y <- datax[!z,]

    Te <- delaunayn(y)
    Ti[[i]] <- Te
    Ts <- tsearchn(y, Te, x)
    indx <- which(!is.na(Ts$idx))
    indnx <- which(is.na(Ts$idx))
    indT <- Ts$idx[indx]
    baryx <- Ts$p[indx,]
    if(!is.matrix(baryx)) baryx <- matrix(baryx,nrow=1)

    for(k in 1:length(p_pcd))
    {
      for(j in unique(indT))
      {
        ind <- which(indT==j)
        xT <- x[indx[ind],]
        yT <- y[Te[j,],]
        baryT <- baryx[ind,]
        if(method=="pe") pcd_simp <- pcd_pe_simp(xT,yT,baryT,r=p_pcd[k])
        if(method=="cs") pcd_simp <- pcd_cs_simp(xT,yT,baryT,tau=p_pcd[k])

        nD    <- indx[ind[pcd_simp$D]]
        D.in[[i,k]] <- append(D.in[[i,k]],nD)
        nSimp  <- pcd_simp$Simp
        Simp[[i,k]] <- c(Simp[[i,k]],nSimp)
      }
    }

    # PCDs, out of convex hull
    if(length(indnx) > 0)
    {
      # find the end triangles and search the points
      ETe <- outer_delaunayn(y)
      indTe <- etsearchn(y,ETe$Te,x[indnx,])

      # PCDs, outside of convex hull
      for(k in 1:length(p_pcd))
      {
        # for all end triangles
        for(j in unique(indTe))
        {
          ind <- which(indTe==j)
          xT <- x[indnx[ind],]
          yT <- ETe$Te[[j]]
          CM <- ETe$CM

          if(method=="pe") pcd.etri <- pcd_pe_Osimp(xT,yT,r=p_pcd[k])
          if(method=="cs") pcd.etri <- pcd_cs_Osimp(xT,yT,CM,tau=p_pcd[k])

          # add new dominating points, their triangle info, their vertex region
          nD    <- indnx[ind[pcd.etri$D]]
          D.out[[i,k]] <- append(D.out[[i,k]],nD)
          nSimp  <- pcd.etri$Simp
          Osimp[[i,k]] <- append(Osimp[[i,k]],nSimp)
        }
      }
    }
  }

  return(list(D.in=D.in,D.out=D.out,C=C,p_pcd=p_pcd,Simp=Simp,Osimp=Osimp,Ti=Ti,classes=cls))
}

#' Classify a set of unlabelled data based on a PCD.
#'
#' @param datan The unlabelled data set
#' @param graph_pcd A proximity catch digraph (PCD)
#' @param p The value of the parameter of the associated proximity map of the PCD
#'
#' @return Predicted labels of the unlabelled data set.
pcd_multiclass_classify <- function(datan,graph_pcd,p)
{
  # distance matrix of points, datan, and dominating sets of pcds
  ncn <- ncol(datan)
  cls <- graph_pcd$classes
  dist_matrix <- NULL
  par <- which(graph_pcd$p.pcd==p)

  # PCD and CCCD, and join them into a single list
  D <- list()
  for(i in 1:length(cls))
    D[[i]] <- c(graph_pcd$D.in[[i,par]],graph_pcd$D.out[[i,par]])
  nc <- unlist(lapply(D,length))
  cl <- rep(cls,nc)

  # full class list
  for(i in 1:length(cls))
  {
    # find pcd triangle distance measures
    Simp <- graph_pcd$Simp[[i,par]]
    if(length(Simp)>0){
      temp <- sapply(Simp,function(simp){
        baryx <- cart2bary(simp,datan)
        if(!is.matrix(baryx)) baryx <- matrix(baryx,nrow=1)
        smallbary <- apply(baryx,1,function(t){
          return(1-(ncn+1)*min(t))
        })
      }, simplify=TRUE)
      dist_matrix <- cbind(dist_matrix,temp)
    }

    # find chcd convex distance measures
    Osimp <- graph_pcd$Osimp[[i,par]]
    if(length(Osimp) > 0){
      temp <- sapply(Osimp,getDisttoOsimp,datan=datan)
      dist_matrix <- cbind(dist_matrix,temp)
    }
  }

  # results
  result <- apply(dist_matrix,1,function(z){
    return(cl[which.min(z)])
  })

  return(result)
}
