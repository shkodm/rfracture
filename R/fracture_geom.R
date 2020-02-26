#' Generate a 3D fracture geometry with trianglar mesh
#' 
#' @param width linear dimension of the square cut of the fracture
#' @param refine refinement level
#' @param power.spectrum power spectrum of the fields (function of frequency)
#' @param corr.profile correlation profile between (function of wave length)
#' @param ... parameters passed to fracture_matrix function
#' 
#' @examples
#' ret = fracture_geom(
#'   width = 1,
#'   refine = 10,
#'   corr.profile = function(lambda) 0.5,
#'   closed=0.1,
#'   power.spectrum = exp_spectrum(scale=0.01, alpha=3)
#' )
#' fracture3d(cut(ret))
#'
#' @export
fracture_geom = function(width=1, refine=1, power.spectrum=exp_spectrum(scale=0.01,alpha=2.5), corr.profile=ogilvie.corr.profile(0.5), ...) {
  n = 6 * refine
  m = 5 * refine
  period = matrix(c(width,0,0,width),2,2)
  ret = fracture_matrix(c(5*n, m), span = matrix(c( 5*width, 0, 3*width, width),2,2), period=period, power.spectrum=power.spectrum, corr.profile=corr.profile, ...)
  
  f1 = ret$f1
  f2 = ret$f2
  A = ret$points[,1]
  dim(A) = ret$dims
  B = ret$points[,2]
  dim(B) = ret$dims
  
  A = rbind(A[(5*n-refine+1):(5*n),1:m]-width*5,     A[1:(n+refine+1),1:m])
  B = rbind(B[(5*n-refine+1):(5*n),1:m]-width*3, B[1:(n+refine+1),1:m])
  f1 = rbind(f1[(5*n-refine+1):(5*n),1:m], f1[1:(n+refine+1),1:m])
  f2 = rbind(f2[(5*n-refine+1):(5*n),1:m], f2[1:(n+refine+1),1:m])
  
  A = cbind(A[,(1*refine+1):(5*refine)],A,A[,1:(refine+1)])
  B = cbind(B[,(1*refine+1):(5*refine)]-width,B,B[,1:(refine+1)]+width)
  f1 = cbind(f1[,(1*refine+1):(5*refine)],f1,f1[,1:(refine+1)])
  f2 = cbind(f2[,(1*refine+1):(5*refine)],f2,f2[,1:(refine+1)])
  
  P = data.frame(x = as.vector(A),y = as.vector(B))
  P$f1 = as.vector(f1)
  P$f2 = as.vector(f2)
  I = 1:length(A)
  dim(I) = dim(A)
  i = rbind(cbind(
    as.vector(I[-nrow(I),-ncol(I)]),
    as.vector(I[-1,-ncol(I)]),
    as.vector(I[-nrow(I),-1])
  ),cbind(
    as.vector(I[-1,-1]),
    as.vector(I[-nrow(I),-1]),
    as.vector(I[-1,-ncol(I)])
  ))
  
  eps = 1e-10
  sel = P$y < width*(1+1/5)+eps & P$y >= width*(-1/5)-eps
  
  ni = rep(0,nrow(P))
  ni[sel] = 1:sum(sel)
  i[] = ni[i]
  i = i[rowSums(i == 0) == 0,]
  P = P[sel,]
  
  P$h = P$f1 - P$f2
  bonds2 = range(P$f1,P$f2)
  cat("Final Bonds:",bonds2[1],bonds2[2],"\n")
  
  ret$width = width
  ret$points=P
  ret$triangles=i
  class(ret) = "fracture_geom"
  return(ret)
}

#' Cut the fracture geometry to a box
#' 
#' @param x fracture_geom object
#' @param eps numerical margin
#' @param ... other arguments
#' @export
cut.fracture_geom = function(x, eps = 1e-9, ...){
  obj = x
  width = obj$width
  snap = function(x) ifelse(x > -eps & x < eps, 0, ifelse(x > width-eps & x < width+eps, width, x))
  i = obj$triangles
  sel = obj$points$x >= -eps & obj$points$x <= width+eps & obj$points$y >= -eps & obj$points$y <= width+eps
  sel = sel[i]
  dim(sel) = dim(i)
  tocut = (! sel[,1]) & sel[,2] & sel[,3]
  obj$points[i[tocut,1],] = (obj$points[i[tocut,1],] + obj$points[i[tocut,3],])/2
  
  sel = obj$points$x >= -eps & obj$points$x <= width+eps & obj$points$y >= -eps & obj$points$y <= width+eps
  ni = rep(0,nrow(obj$points))
  ni[sel] = 1:sum(sel)
  i[] = ni[i]
  i = i[rowSums(i == 0) == 0,]
  ret = obj
  ret$points = obj$points[sel,]
  ret$points$x = snap(ret$points$x)
  ret$points$y = snap(ret$points$y)
  ret$triangles = i
  ret
}

#' Calculate the volume of fracture geometry
#' 
#' @param x the fracture_geom object to calculate volume of
#' @param ... other parameters (ignored)
#' 
#' @export
volume.fracture_geom = function(x, ...) {
  obj = x
  v = obj$points[,c("x","y")]
  v1 = v[obj$triangles[,2],] - v[obj$triangles[,1],]
  v2 = v[obj$triangles[,3],] - v[obj$triangles[,1],]
  a = abs(1/2*(v1[,1]*v2[,2] - v1[,2]*v2[,1]))
  h = 1/3*(obj$points$h[obj$triangles[,1]] + obj$points$h[obj$triangles[,2]] + obj$points$h[obj$triangles[,3]])
  sum(a*h)
}

#' Extract the overlapping part of two surfaces of a fracture
#' 
#' @param obj the fracture_geom object
#' @param touch what to do with the touching part (see description)
#' @description 
#' The possible values of touch are: "exclude" - exclude touching part, "include" - leave untouched, "only" - select only the touching part.
#' 
#' @export
touching = function(obj,touch="exclude") {
  P = obj$points
  i = obj$triangles
  sel = P$f1 == P$f2
  
  sel = sel[i]
  dim(sel) = dim(i)
  sel = rowSums(sel) != 3
  if (touch == "exclude") {
    
  } else if (touch == "include") {
    sel[] = TRUE
  } else if (touch == "only") {
    sel = ! sel
  } else stop("invalid value for touch")
  i = i[sel,]
  sel = rep(FALSE,nrow(P))
  sel[i[,1]] = TRUE
  sel[i[,2]] = TRUE
  sel[i[,3]] = TRUE
  
  ni = rep(0,nrow(P))
  ni[sel] = 1:sum(sel)
  i[] = ni[i]
  i = i[rowSums(i == 0) == 0,]
  P = P[sel,]
  obj$points=P
  obj$triangles=i
  return(obj)
}

