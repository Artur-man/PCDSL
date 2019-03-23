#' Find the pcd of a single simplex for a given parameter r.
#'
#' @param x The set of points of the target class
#' @param y The coordinates of the vertices of the n-simplex
#' @param bary The barycentric coordinates of the points of x with respect to the n-simplex
#' @param tau The CS-PCD parameter
#'
#' @return The dominating set of the n-simplex and the associated proximity regions
pcd_cs_simp <- function(x,y,bary,tau)
{
  if(!is.matrix(bary)) bary <- matrix(bary,nrow=1)
  if(!is.matrix(x)) x <- matrix(x,nrow=1)
  nx <- nrow(x)
  macheps <- .Machine$double.eps
  allind <- 1:nrow(y)

  ind <- apply(bary,1,which.min)
  uniqreg <- sort(unique(ind))
  med <- apply(y,2,mean)
  Simp <- list()

  disty <- rep(0,nrow(y))
  for(j in allind){
    disty[j] <- dist_to_plane(med,y[-j,])
  }

  for(i in 1:nx)
  {
    simpx <- NULL
    distx <- dist_to_plane(x[i,],y[-ind[i],])
    simx <- tau*(distx/disty[ind[i]])
    for(j in allind){
      line  <- get_line_par_to(x[i,],rbind(med,y[j,]))
      simpx <- rbind(simpx,x[i,]+(line[2,]-line[1,])*simx)
    }
    Simp <- c(Simp,list(simpx))
  }

  if(tau > 1)
  {
    for(k in 1:length(Simp)){
      simp <- Simp[[k]]
      for(i in allind){
        dist1 <- dist_to_plane(simp[i,],y[-i,])
        dist2 <- dist_to_plane(simp[i,],simp[-i,])
        if(dist2 > dist1) {
          sim <- dist1/dist2
          temp <- apply(simp[-i,],1,function(t){
            return(simp[i,]+(t-simp[i,])*sim)
          })
          simp[-i,] <- t(temp)
        }
      }
      Simp[[k]] <- simp
    }
  }

  A <- in_simp(Simp,x)
  A <- t(A)
  D <- dominate_greedy_matrix(A)
  Dsimp <- Simp[D]

  return(list(D=D,Simp=Dsimp))
}

#' Find the pcd of an outer simplex for a given parameter r.
#'
#' @param x The set of points of the target class
#' @param y The coordinates of the vertices of the outer simplex
#' @param CM The median point of the convex hull of the non-target class points
#' @param tau The CS-PCD parameter
#'
#' @return The dominating point of the outer simplex and the associated proximity region
pcd_cs_Osimp <- function(x,y,CM,tau)
{
  if(!is.matrix(x)) x <- matrix(x,nrow=1)

  nx <- nrow(x)
  nc <- ncol(y)
  macheps <- .Machine$double.eps
  allind <- 1:nc
  Simp <- list()

  c_y <- t(combn(allind,nc-1))
  c_y <- cbind(c_y,c_y[,1]+nc)
  c_y <- rbind(allind,c_y)
  distx <- apply(c_y,1,function(t){
                       dist_to_plane(x,y[t,])
  })
  distx <- matrix(distx,ncol=nc+1)
  ind <- apply(distx,1,which.min)

  inc <- get_incenter_Osimp(y[1:nc,],CM)
  disty <- dist_to_plane(inc,y[allind,])

  inc_f <- get_plane_par_to(inc,y[allind,])
  y_n <- y
  for(j in allind){
    y_line <- rbind(y[j,],y[j+nc,])
    inc_p <- get_line_plane_inter(y_line,inc_f)
    distp <- sqrt(sum((inc_p-y[j,])^2))
    y_n[j+nc,] <- y[j,] + (inc_p-y[j,])*2
  }

  Simp <- list()
  for(i in 1:nx)
  {
    simpx <- matrix(nrow=nc*2,ncol=nc)
    simx <- tau*(distx[i,ind[i]]/disty)
    for(j in allind){
      line  <- get_line_par_to(x[i,],rbind(inc,y_n[j,]))
      simpx[j,] <- x[i,]+(line[2,]-line[1,])*simx
      line  <- get_line_par_to(x[i,],rbind(inc,y_n[j+nc,]))
      simpx[j+nc,] <- x[i,]+(line[2,]-line[1,])*simx
    }
    Simp <- c(Simp,list(simpx))
  }

  if(tau > 1)
  {
    for(k in 1:length(Simp)){
      simp <- Simp[[k]]
      for(i in 1:nrow(c_y)){
        if(all(c_y[i,]==allind)){
          c_y_c <- allind+nc
          dist1 <- dist_to_plane(simp[nc+1,],y[allind,])
          dist2 <- dist_to_plane(simp[nc+1,],simp[allind,])
          if(dist2 > dist1) {
            sim <- as.double(dist1/dist2)
            simp[allind,] <- simp[c_y_c,]+(simp[allind,]-simp[c_y_c,])*sim
          }
        } else {
          c_y_c <- setdiff(allind,c_y[i,-nc])
          c_y_c <- c(c_y_c,c_y_c+nc)
          dist1 <- dist_to_plane(x[k,],y[c_y[i,],])
          dist2 <- dist_to_plane(x[k,],simp[c_y[i,],])
          if(dist2 > dist1){
            # one side first
            dist1 <- dist_to_plane(simp[c_y_c[1],],y[c_y[i,],])
            dist2 <- dist_to_plane(simp[c_y_c[1],],simp[c_y[i,],])
            sim <- as.double(dist1/dist2)
            c_y_c_f_1 <- setdiff(allind,c_y_c[1])
            temp1 <- matrix(simp[c_y_c_f_1,],ncol=nc)
            temp1 <- apply(temp1,1,function(t){
              return(simp[c_y_c[1],]+(t-simp[c_y_c[1],])*sim)
            })
            temp1 <- t(temp1)

            dist1 <- dist_to_plane(simp[c_y_c[2],],y[c_y[i,],])
            dist2 <- dist_to_plane(simp[c_y_c[2],],simp[c_y[i,],])
            sim <- as.double(dist1/dist2)
            c_y_c_f_2 <- c_y_c_f_1+nc
            temp2 <- matrix(simp[c_y_c_f_2,],ncol=nc)
            temp2 <- apply(temp2,1,function(t){
              return(simp[c_y_c[2],]+(t-simp[c_y_c[2],])*sim)
            })
            temp2 <- t(temp2)

            simp[c(c_y_c_f_1,c_y_c_f_2),] <- rbind(temp1,temp2)
          }
        }
      }
      Simp[[k]] <- simp
    }
  }

  # adjacency matrix of the CS-PCD and its dominating set
  A <- in_Osimp(Simp,x)
  A <- t(A)
  D <- dominate_greedy_matrix(A)
  Dsimp <- Simp[D]

  return(list(Simp=Dsimp,D=D))
}
