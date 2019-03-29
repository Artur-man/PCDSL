#' Calculate the distance from a set of points to a plane
#'
#' @param p m-by-d matrix in which each row is the Cartesian coordinates of a point.
#' @param plane d-by-d matrix representing a plane.
#'
#' @return The set of distances between each point in p and a plane
dist_to_plane <- function(p,plane)
{
  n <- get_normal_vector(plane)

  if(is.matrix(p)) {
    distx <- apply(p,1,function(x){
      tmp <- x-plane[1,]
      return(abs(crossprod(tmp,n)))
    })
  } else {
    tmp <- p-plane[1,]
    distx <- abs(crossprod(tmp,n))
  }

  return(distx)
}

#' Find if points are in a set of d-simplices
#'
#' @param simp A list of d+1-by-d matrices each representing a d-simplex.
#' @param x The Cartesian coordinates of a reference point.
#'
#' @return A logical vector indicating in which simplex the reference point x resides.
in_simp <- function(simp,x)
{
  if(!is.matrix(x)) x <- matrix(x,ncol=length(x))
  which_simp <- sapply(simp,function(z){
    baryx <- cart2bary(z,x)
    in_simp_x <- apply(baryx,1,function(t){
      if(any(t < 0)){
        return(FALSE)
      } else return(TRUE)
    })
    return(in_simp_x)
  },simplify = TRUE)
  return(which_simp)
}

#' Find if points are in a set of outer simplices
#'
#' @param Osimp A list of d+1-by-d matrices each representing an outer simplex.
#' @param x The Cartesian coordinates of a reference point.
#'
#' @return A logical vector indicating in which simplex the reference point x resides.
in_Osimp <- function(Osimp,x)
{
  if(!is.matrix(x)) x <- matrix(x,ncol=length(x))
  which_simp <- sapply(Osimp,function(z){
    baryx <- cart2genbary(z,x)
    in_Osimp_x <- apply(baryx,1,function(t){
      if(all(t > 0 & t < 1)){
        return(TRUE)
      } else return(FALSE)
    })
    return(in_Osimp_x)
  },simplify = TRUE)
  return(which_simp)
}

#' Find a parallel plane passing a reference point.
#'
#' @param x The Cartesian coordinates of a reference point
#' @param plane d-by-d matrix representing a plane.
#'
#' @return A plane passing a reference point x which is parallel to the plane.
get_plane_par_to <- function(x,plane){

  A <- plane[-1,]
  if(!is.matrix(A)) A <- matrix(A,nrow=1)
  A <- apply(A,1,function(t){return(t-plane[1,])})
  P <- A%*%(solve(t(A)%*%A))%*%t(A)
  yx <- matrix(x-plane[1,],ncol=1)
  newx <- plane[1,] + P%*%yx
  newy <- t(apply(plane,1,function(t){return((t-t(newx)+x))}))
  return(newy)
}

#' Find a parallel ray passing a reference point.
#'
#' @param x The Cartesian coordinates of a reference point
#' @param ray 2-by-d matrix representing a ray.
#'
#' @return A plane passing a reference point x which is parallel to the ray.
get_line_par_to <- function(x,ray)
{
  vecray <- ray[2,]-ray[1,]
  return(rbind(x,x+vecray))
}

#' Find the point of intersection between a ray and a plane.
#'
#' @param ray 2-by-d matrix representing a ray.
#' @param plane d-by-d matrix representing a plane.
#'
#' @return A point on a ray.
get_line_plane_inter <- function(ray,plane)
{
  temp <- plane[-1,]
  if(!is.matrix(temp)) temp <- matrix(temp,nrow=1)
  line <- matrix(ray[2,]-ray[1,],ncol=1)
  A <- cbind(line,apply(temp,1,function(t){return(t-plane[1,])}))
  b <- plane[1,]-ray[1,]
  a <- solve(A,b)
  return((1-a[1])*ray[1,]+a[1]*ray[2,])
}


#' Get the volume of a d-simplex.
#'
#' @param simp d+1-by-d matrix representing a simplex.
#'
#' @return The volume of the simplex.
get_plane_volume <- function(simp){

  vecs <- t(apply(simp,1,function(x){
    return(x-simp[1,])
  }))
  n.simp <- get_normal_vector(vecs)
  vecs <- rbind(n.simp,vecs)
  dummy <- rep(1,nrow(vecs))
  vecs <- cbind(dummy,vecs)

  return(abs(det(vecs)))
}

#' Get the incenter of an outer simplex
#'
#' @param facet d-by-d matrix representing a facet of the convex hull of the non-target class
#' @param m The median point of the convex hull of the non-target class
#'
#' @return The incenter of the outer simplex
get_incenter_Osimp <- function(facet,m){

  cm <- get_incenter_simp(rbind(m,facet))
  distcm <- dist_to_plane(cm,facet)
  distm <- dist_to_plane(m,facet)
  distx <- distm/(distm/distcm-2)
  distt <- (distm+distx)/(distm-distcm)
  distt <- as.double(distt)
  Im <- (cm-m)*distt  + m

  return(Im)
}

#' Get the incenter of an n-simplex
#'
#' @param simp d+1-by-d matrix representing an d-simplex
#'
#' @return The incenter of the outer simplex
get_incenter_simp <- function(simp){

  nc <- ncol(simp)
  allind <- 1:nrow(simp)
  c_ind <- t(combn(allind,nc))
  area_simp <- apply(c_ind,1,function(t){
                            return(get_plane_volume(simp[t,]))
  })
  area_simp <- rev(area_simp)
  bary_simp <- area_simp/sum(area_simp)
  bary_simp <- matrix(bary_simp,nrow=1)

  cm <- bary2cart(simp,bary_simp)

  return(cm)
}

#' Get the normal vector of a plane
#'
#' @param plane d-by-d matrix representing a plane
#'
#' @return The normal vector of the plane.
get_normal_vector <- function(plane){

  if(!is.matrix(plane)) plane <- matrix(plane,nrow=1)
  if(nrow(plane) > ncol(plane)) plane <- plane[1:ncol(plane),]
  veclen <- nrow(plane)-1
  vecpla <- plane[-1,]-outer(rep(1,veclen),plane[1,])
  A <- cbind(t(vecpla),runif(ncol(plane)))
  Q <- qr.Q(qr(A))
  n <- Q[,veclen+1]
  n <- n/sqrt(sum(n^2))

  return(n)
}

#' Find the index of a list of outer simplices that contains a set of points x.
#'
#' @param y The set of points of the non-target class
#' @param Osimp The list of the outer simplices of the convex hull of the non-target class.
#' @param x The set of reference points of the target class
#'
#' @return A vector indicating the indices of the list of outer simplices for each point of x.
etsearchn <- function(y,Osimp,x){

  indt <- (ncol(y)+1):(2*ncol(y))

  tflag <- NULL
  for(i in 1:length(Osimp)){
    cartx <- cart2genbary(Osimp[[i]],x)
    flag <- apply(cartx,1,function(t){
                         return(all(t[indt]>0))
    })
    tflag <- cbind(tflag,flag)
  }

  indices <- apply(tflag,1,function(t){
                          if(any(t)) return(which(t))
                          else return(F)
  })

  return(indices)
}

#' The list of the outer simplices of the Delaunay tessellation of the non-target class
#'
#' @param y The set of points of the non-target class
#'
#' @return The list of outer simplices of the Delaunay tessellation
outer_delaunayn <- function(y){

  cv <- convhulln(y)
  bound.vec <- unique(as.vector(cv))
  m <- apply(y[bound.vec,],2,mean)
  CM <- m

  rays <- apply(y[bound.vec,],1,function(t){
    tp <- t-m
    tp <- tp/sqrt(sum(tp^2))
    return(tp)
  })
  rays <- t(rays)

  yt <- NULL
  for(i in 1:nrow(rays)){
    yt <- rbind(yt,y[bound.vec[i],]+rays[i,])
  }

  Osimp <- list()
  for(i in 1:nrow(cv)){
    m <- match(cv[i,],bound.vec)
    Osimp[[i]] <- rbind(y[cv[i,],],yt[m,])
  }
  return(list(Te=Osimp,CM=CM))
}

#' Given the Cartesian coordinates of a set of points,
#' compute the generalized barycentric coordinates of these
#' points with respect to an outer simplex using the method of Warren 1996.
#'
#' @param Osimp A 2d-by-d matrix representing a outer simplex
#' @param x A set of points of the target class
#'
#' @return The barycentric coordinates of the each points of x.
#'
#' @references
#'
#' Joe Warren. Barycentric Coordinates for Convex Polytopes. Advances in Computational Mathematics, 6:97–108, 1996.
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
  face_list <- apply(faces,2,list)
  list_face <- c(list_face,lapply(face_list,unlist))
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

#' Get the distances between a set of points and a outer simplex.
#'
#' @param Osimp A 2d-by-d matrix representing a outer simplex.
#' @param x The set of points of the target class
#'
#' @return The set of distances between each point of x and a outer simplex
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



