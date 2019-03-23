#' Calculate the distance from a set of points to a plane
#'
#' @param p M-by-N matrix in which each row is the Cartesian coordinates of a point.
#' @param pla N-by-N matrix representing a plane.
#'
#' @return The set of distances between each point in p and plane pla.
dist_to_plane <- function(p,pla)
{
  n <- get_normal_vector(pla)

  if(is.matrix(p)) {
    distx <- apply(p,1,function(x){
      tmp <- x-pla[1,]
      return(abs(crossprod(tmp,n)))
    })
  } else {
    tmp <- p-pla[1,]
    distx <- abs(crossprod(tmp,n))
  }

  return(distx)
}

#' find if points are in a set of n-simplexes
#'
#' @param simp A list of N+1-by-N matrices each representing an n-simplex.
#' @param x The Cartesian coordinates of a reference point
#'
#' @return A logical vector indicating in which simplex the reference point x resides.
in_simp <- function(simp,x)
{
  if(!is.matrix(x)) x <- matrix(x,ncol=length(x))
  matx <- sapply(simp,function(z){
    baryx <- cart2bary(z,x)
    smallbary <- apply(baryx,1,function(t){
      if(any(t < 0)){
        return(FALSE)
      } else return(TRUE)
    })
    return(smallbary)
  },simplify = TRUE)
  return(matx)
}

#' find if points are in a set of outer simplexes
#'
#' @param Osimp A list of N+1-by-N matrices each representing an outer simplex.
#' @param x The Cartesian coordinates of a reference point
#'
#' @return A logical vector indicating in which simplex the reference point x resides.
in_Osimp <- function(Osimp,x)
{
  if(!is.matrix(x)) x <- matrix(x,ncol=length(x))
  matx <- sapply(Osimp,function(z){
    baryx <- cart2genbary(z,x)
    smallbary <- apply(baryx,1,function(t){
      if(all(t > 0 & t < 1)){
        return(TRUE)
      } else return(FALSE)
    })
    return(smallbary)
  },simplify = TRUE)
  return(matx)
}

#' Find a parallel plane passing a point, by means of orthogonal projection
#'
#' @param x The Cartesian coordinates of a reference point
#' @param pla N-by-N matrix representing a plane.
#'
#' @return A plane passing the point x which the parallel to the plane pla.
get_plane_par_to <- function(x,pla){

  A <- pla[-1,]
  if(!is.matrix(A)) A <- matrix(A,nrow=1)
  A <- apply(A,1,function(t){return(t-pla[1,])})
  P <- A%*%(solve(t(A)%*%A))%*%t(A)
  yx <- matrix(x-pla[1,],ncol=1)
  newx <- pla[1,] + P%*%yx
  newy <- t(apply(pla,1,function(t){return((t-t(newx)+x))}))
  return(newy)
}

#' Find a parallel ray passing the point, by means of orthogonal projection
#'
#' @param x The Cartesian coordinates of a reference point
#' @param ray 2-by-N matrix representing a ray.
#'
#' @return A plane passing the point x which the parallel to the plane pla.
get_line_par_to <- function(x,ray)
{
  vecray <- ray[2,]-ray[1,]
  return(rbind(x,x+vecray))
}

# find the point where the line and the plane intersects

#' Find the point where the ray and the plane intersects
#'
#' @param ray 2-by-N matrix representing a ray.
#' @param pla N-by-N matrix representing a plane.
#'
#' @return A plane passing the point x which the parallel to the plane pla.
get_line_plane_inter <- function(ray,pla)
{
  temp <- pla[-1,]
  if(!is.matrix(temp)) temp <- matrix(temp,nrow=1)
  line <- matrix(ray[2,]-ray[1,],ncol=1)
  A <- cbind(line,apply(temp,1,function(t){return(t-pla[1,])}))
  b <- pla[1,]-ray[1,]
  a <- solve(A,b)
  return((1-a[1])*ray[1,]+a[1]*ray[2,])
}


#' Get the volume of the simplex.
#'
#' @param simp N-by-N-1 matrix representing a simplex.
#'
#' @return The volume of the simplex.
get_plane_volume <- function(simp){

  mat <- t(apply(simp,1,function(x){
    return(x-simp[1,])
  }))
  n.simp <- get_normal_vector(mat)
  mat <- rbind(n.simp,mat)
  temp <- rep(1,nrow(mat))
  mat <- cbind(temp,mat)

  return(abs(det(mat)))
}

#' Get the incenter of an outer simplex
#'
#' @param pla N-by-N matrix representing a facet of the convex hull of the non-target class
#' @param m The median point of the convex hull of the non-target class
#'
#' @return The incenter of the outer simplex
get_incenter_Osimp <- function(pla,m){

  cm <- get_incenter_simp(rbind(m,pla))
  distcm <- dist_to_plane(cm,pla)
  distm <- dist_to_plane(m,pla)
  distx <- distm/(distm/distcm-2)
  temp <- (distm+distx)/(distm-distcm)
  temp <- as.double(temp)
  Im <- (cm-m)*temp  + m

  return(Im)
}

#' Get the incenter of an n-simplex
#'
#' @param simp N+1-by-N matrix representing an n-simplex
#'
#' @return The incenter of the outer simplex
get_incenter_simp <- function(simp){

  nc <- ncol(simp)
  allind <- 1:nrow(simp)
  c.ind <- t(combn(allind,nc))
  area.simp <- apply(c.ind,1,function(t){
                            return(get_plane_volume(simp[t,]))
  })
  # the bary. coord are reversed cos of c.ind
  area.simp <- rev(area.simp)
  bary.simp <- area.simp/sum(area.simp)
  bary.simp <- matrix(bary.simp,nrow=1)

  cm <- bary2cart(simp,bary.simp)

  return(cm)
}

#' Get the normal vector of a plane
#'
#' @param pla N-by-N matrix representing a plane
#'
#' @return The normal vector of the plane.
get_normal_vector <- function(pla){

  if(!is.matrix(pla)) pla <- matrix(pla,nrow=1)

  # check if more than necessary points are presented
  # if so, erase some and keep enough of them
  if(nrow(pla) > ncol(pla)) pla <- pla[1:ncol(pla),]

  # vectors of non-target points in the list in terms of points p0
  veclen <- nrow(pla)-1
  vecpla <- pla[-1,]-outer(rep(1,veclen),pla[1,])

    # find the normal vector by solving Yn=0,
  # where Y is the vector and n is the normal
  A <- cbind(t(vecpla),runif(ncol(pla)))
  Q <- qr.Q(qr(A))
  n <- Q[,veclen+1]

  # normalize
  n <- n/sqrt(sum(n^2))

  return(n)
}

#' Search for the outer simplicies that a set of points belong to
#'
#' @param y The set of points of the non-target class
#' @param Osimp The indices of the facets associated with the outer simplicies
#' @param x The set of reference points of the target class
#'
#' @return A vector indicating the indices of Osimp in which each point of x is found
etsearchn <- function(y,Osimp,x){

  # if the barycentric coordinates of the matching vertices are satisfied
  # then the point is in that region
  indt <- (ncol(y)+1):(2*ncol(y))

  tflag <- NULL
  for(i in 1:length(Osimp)){
    cartx <- cart2genbary(Osimp[[i]],x)
    flag <- apply(cartx,1,function(t){
                         return(all(t[indt]>0))
    })
    tflag <- cbind(tflag,flag)
  }

  result <- apply(tflag,1,function(t){
                          if(any(t)) return(which(t))
                          else return(F)
  })

  return(result)
}


#' The list of the outer simplicies of the Delaunay tessellation of the non-target class
#'
#' @param y The set of points of the non-target class
#'
#' @return The list of outer simplicies of the Delaunay tessellation
outer_delaunayn <- function(y){

  # find the boundary vertices and the mean of the convex hull
  cv <- convhulln(y)
  bound.vec <- unique(as.vector(cv))
  m <- apply(y[bound.vec,],2,mean)
  CM <- m

  # normalize and then tak the transpose
  rays <- apply(y[bound.vec,],1,function(t){
    tp <- t-m
    tp <- tp/sqrt(sum(tp^2))
    return(tp)
  })
  rays <- t(rays)

  # get end points of rays
  yt <- NULL
  for(i in 1:nrow(rays)){
    yt <- rbind(yt,y[bound.vec[i],]+rays[i,])
  }

  # get the end triangles
  Osimp <- list()
  for(i in 1:nrow(cv)){
    m <- match(cv[i,],bound.vec)
    Osimp[[i]] <- rbind(y[cv[i,],],yt[m,])
  }
  return(list(Te=Osimp,CM=CM))
}

#' Given the Cartesian coordinates of one or more points, compute the generalized barycentric coordinates of these
#' points with respect to an outer simplex using the method of Warren 1996.
#'
#' @param Osimp the outer
#' @param x A set of points of the target class
#'
#' @return The barycentric coordinates of the each points of x.
cart2genbary <- function(Osimp,x){

  if(!is.matrix(x)) x <- matrix(x,nrow=1)
  bound_vec <- 1:nrow(Osimp)
  o_face <- (ncol(Osimp)+1):(2*ncol(Osimp))
  face1 <- combn(1:ncol(Osimp),ncol(Osimp)-1)
  face2 <- apply(face1,2,function(t){
                         return(o_face[t])
  })
  list_face <- list(1:ncol(Osimp),(ncol(Osimp)+1):(2*ncol(Osimp)))
  faces <- rbind(face1,face2)
  temp <- apply(faces,2,list)
  list_face <- c(list_face,lapply(temp,unlist))
  norm_vecs <- sapply(list_face,function(t){
    n <- get_normal_vector(Osimp[t,])
    z <- setdiff(bound_vec,t)
    if((Osimp[t[2],]-Osimp[z[1],]) %*% n < 0) n <- -n
    return(n)
  },simplify = TRUE)
  norm_vecs <- t(norm_vecs)
  baryx <- matrix(0,nrow=nrow(x),ncol=length(bound_vec))
  for(i in bound_vec){
    flag <- sapply(list_face,function(t){
                       return(any(t==i))
    },simplify = TRUE)
    ind_n <- which(flag)
    det_n <- abs(det(norm_vecs[ind_n,]))
    res <- apply(norm_vecs[ind_n,],1,function(t){
                                     prod_n <- apply(x,1,function(z){
                                                         return(t%*%(Osimp[i,]-z))
                                     })
                                     return(prod_n)
    })
    if(!is.matrix(res)) res <- matrix(res,nrow=1)
    denom_n <- apply(res,1,prod)
    baryx[,i] <- det_n/denom_n
  }
  baryx <- t(apply(baryx,1,function(t){
                         return(t/sum(t))
  }))

  return(baryx)
}

#' Get the distances of a list of points from an outer simplex.
#'
#' @param Osimp A 2N-by-N matrix representing an outer simplex.
#' @param x The set of points of the target class
#'
#' @return The set of distances for each point of x
getDisttoOsimp <- function(Osimp,x){

  nc <- ncol(Osimp)
  ind <- 1:nc
  m <- apply(Osimp,2,mean)

  cb <- t(combn(nc,nc-1))
  cb <- cbind(cb,cb[,1]+nc)
  ch <- rbind(ind,ind+nc,cb)

  alldist <-  apply(x,1,
                    function(z){
                      yn <- rep(0,nrow(ch))
                      for(i in 1:nrow(ch))
                      {
                        chdata <- Osimp[ch[i,],]
                        temp <- chdata[1,]
                        chdata <- apply(chdata,1,function(t){
                          return(chdata[1,]-t)
                        })
                        chdata <- chdata[,-1]
                        A <- cbind(z-m,chdata)
                        b <- temp-m
                        y <- solve(A,b)
                        yn[i] <- y[1]
                      }
                      if(min(yn[yn>=0])==Inf) print(yn)
                      res <- 1/min(yn[yn>=0])
                      return(res)
                    })
  return(alldist)
}



