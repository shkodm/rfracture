#' Generate a 3D fracture geometry with trianglar mesh
#' 
#' @param width linear dimension of the square cut of the fracture
#' @param refine refinement level
#' @param ... parameters passed to fracture_matrix function
#' 
#' @examples
#' ret = fracture_geom(
#'   width = 1,
#'   refine = 50,
#'   corr.profile = function(lambda) 0.9,
#'   gap = 0.1,
#'   power.spectrum = exp_spectrum(scale=0.004, fractal.dimension=3-0.5)
#' )
#' fracture3d(ret)
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
#' @param obj fracture_geom object
#' @param eps numerical margin
#' 
#' @export
cut.fracture_geom = function(obj, eps = 1e-9){
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

