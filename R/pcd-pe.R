#' Find the pcd of a single simplex for a given parameter r.
#'
#' @param x The set of points of the target class  
#' @param y The coordinates of the vertices of the n-simplex
#' @param bary The barycentric coordinates of the points of x with respect to the n-simplex 
#' @param r The PE-PCD parameter 
#'
#' @return The dominating set of the n-simplex and the associated proximity regions
pcd_pe_simp <- function(x,y,bary,r)
{
  if(!is.matrix(bary)) bary <- matrix(bary,nrow=1) 
  if(!is.matrix(x)) x <- matrix(x,nrow=1) 

  ind <- apply(bary,1,which.max) 
  nx <- nrow(x)                  
  macheps <- .Machine$double.eps 
  uniqreg <- sort(unique(ind))   

  extlist <- NULL   
  Tri <- list()
  A <- NULL 
  
  for(i in uniqreg)
  {
    indx <- which(ind==i)              
    disty <- dist_to_plane(y[i,],y[-i,]) 
    distxall <- disty - dist_to_plane(x,y[-i,])  
    
    dom <- which.max(distxall[indx])
    distx <- distxall[indx[dom]]         
    distt <- r*distx                    
    if(distt < distx) distt=distx        
    if(distt > disty) distt=disty        
    
    A <- rbind(A,distxall <= as.vector(distt)) 
    extlist <- append(extlist,indx[dom]) 
    ratio <- as.vector(distt/disty)
    Prox <- t(apply(y,1,function(z){
      return(y[i,] + ratio*(z-y[i,]))
    }))
    Tri <- c(Tri,list(Prox))
  } 
  
  powdom <- powerset(1:length(uniqreg))     
  ordpowdom <- order(sapply(powdom,length)) 
  powdom <- powdom[ordpowdom]
  
  for(i in powdom)
  {
    adj <- A[i,] 
    
    adj <- matrix(adj,nrow=length(i))
    alladj <- apply(adj,2,any)
    
    if(all(alladj==TRUE)) {
      D <- extlist[i]
      DTri <- Tri[i]
      break
    }
  }
    
  return(list(D=D,Tri=DTri))
}

#' Find the pcd of an outer simplex for a given parameter r.
#'
#' @param x The set of points of the target class  
#' @param y The coordinates of the vertices of the outer simplex 
#' @param r The PE-PCD parameter 
#'
#' @return The dominating point of the outer simplex and the associated proximity region
pcd_pe_Osimp <- function(x,y,r)
{
  if(!is.matrix(x)) x <- matrix(x,nrow=1) 
  
  nc <- ncol(y)
  face <- y[1:nc,]
  o.face <- y[(nc+1):(2*nc),]
  distx <- apply(x,1,dist_to_plane,pla=face)
  disty <- apply(o.face,1,dist_to_plane,pla=face)
  dom <- which.max(distx)
  
  if(r<1) r=1 
  ray <- o.face - face
  ratio <- (r*distx[dom])/disty
  oppos.face <- face + ratio*ray
  
  D <- dom
  DTri <- list(rbind(face,oppos.face))
  return(list(D=D,Tri=DTri))
}