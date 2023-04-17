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
#'   power.iso = exp_spectrum(scale=0.01, alpha=3)
#' )
#' if (require(rgl)) fracture3d(cut(ret))
#'
#' @export
fracture_geom = function(width=1, refine=1, method=c("triangles","diagonals"), period = width, ...) {
  method = match.arg(method)
  if (length(period) == 1) period = rep(period,2)
  if (length(period) == 2) period = matrix(c(period[1],0,0,period[2]),2,2)
  if (!(inherits(period, "matrix") && identical(dim(period),c(2L,2L)))) stop("Wrong period argument")
  
  if (method == "triangles") {
    n = 6 * refine
    m = 5 * refine
    dims = c(5*n, m)
    span = matrix(c( 5*width, 0, 3*width, width),2,2)
  } else if (method == "diagonals") {
    n = 5 * refine
    m = 5 * refine
    dims = c(n, m)
    span = matrix(c(width, 0, 0, width),2,2)
  }
  ret = fracture_matrix(dims, span = span, period=period, ...)
  
  f1 = ret$f1
  f2 = ret$f2
  A = ret$points[,1]
  dim(A) = ret$dims
  B = ret$points[,2]
  dim(B) = ret$dims
  
  if (method == "triangles") {
    A = rbind(A[(5*n-refine+1):(5*n),1:m]-width*5,     A[1:(n+refine+1),1:m])
    B = rbind(B[(5*n-refine+1):(5*n),1:m]-width*3, B[1:(n+refine+1),1:m])
    f1 = rbind(f1[(5*n-refine+1):(5*n),1:m], f1[1:(n+refine+1),1:m])
    f2 = rbind(f2[(5*n-refine+1):(5*n),1:m], f2[1:(n+refine+1),1:m])
    
    A = cbind(A[,(1*refine+1):(5*refine)],A,A[,1:(refine+1)])
    B = cbind(B[,(1*refine+1):(5*refine)]-width,B,B[,1:(refine+1)]+width)
    f1 = cbind(f1[,(1*refine+1):(5*refine)],f1,f1[,1:(refine+1)])
    f2 = cbind(f2[,(1*refine+1):(5*refine)],f2,f2[,1:(refine+1)])
  } else if (method == "diagonals") {
    A  = rbind(A[(n-refine+1):(n),]-width, A, A[1:(refine+1),]+width)
    B  = rbind( B[(n-refine+1):(n),],  B,  B[1:(refine+1),])
    f1 = rbind(f1[(n-refine+1):(n),], f1, f1[1:(refine+1),])
    f2 = rbind(f2[(n-refine+1):(n),], f2, f2[1:(refine+1),])
    
    A  = cbind( A[,(n-refine+1):(n)],  A,  A[,1:(refine+1)])
    B  = cbind( B[,(n-refine+1):(n)]-width, B, B[,1:(refine+1)]+width)
    f1 = cbind(f1[,(n-refine+1):(n)], f1, f1[,1:(refine+1)])
    f2 = cbind(f2[,(n-refine+1):(n)], f2, f2[,1:(refine+1)])
  } else stop("unknown method")
  P = data.frame(x = as.vector(A),y = as.vector(B))
  P$f1 = as.vector(f1)
  P$f2 = as.vector(f2)
  I = 1:length(A)
  dim(I) = dim(A)
  if (method == "triangles") {
    sel = rep(TRUE, (nrow(I) - 1) * (ncol(I) - 1))
  } else if (method == "diagonals") {
    sel = sample(c(TRUE,FALSE),size = (nrow(I) - 1) * (ncol(I) - 1), replace = TRUE)
  }
  i = rbind(cbind(
    as.vector(I[-nrow(I),-ncol(I)]),
    as.vector(I[-1,-ncol(I)]),
    as.vector(I[-nrow(I),-1])
  )[sel,],cbind(
    as.vector(I[-1,-1]),
    as.vector(I[-nrow(I),-1]),
    as.vector(I[-1,-ncol(I)])
  )[sel,],cbind(
    as.vector(I[-1,-1]),
    as.vector(I[-nrow(I),-1]),
    as.vector(I[-nrow(I),-ncol(I)])
  )[!sel,],cbind(
    as.vector(I[-nrow(I),-ncol(I)]),
    as.vector(I[-1,-ncol(I)]),
    as.vector(I[-1,-1])
  )[!sel,])
  
  eps = 1e-10
  sel = P$y < width*(1+1/5)+eps & P$y >= width*(-1/5)-eps
  
  ni = rep(0,nrow(P))
  ni[sel] = 1:sum(sel)
  i[] = ni[i]
  i = i[rowSums(i == 0) == 0,]
  P = P[sel,]
  
  P$h = P$f1 - P$f2
  P$fm = (P$f1 + P$f2)/2
  bonds2 = range(P$f1,P$f2)
#  cat("Final Bonds:",bonds2[1],bonds2[2],"\n")
  
  ret$f1 = NULL
  ret$f2 = NULL
  ret$width = width
  ret$points=P
  ret$triangles=i
  ret$edge = matrix(nrow=0,ncol=2)
  ret$vertex = matrix(nrow=0,ncol=1)
  class(ret) = "fracture_geom"
  return(ret)
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

#' Summary of surfaces
#' 
#' @param obj fracture_geom object
#' @param surface the variables over which to gather statistics
#' 
#' @export
surface_summary = function(obj, surface=c("f1","f2","h","fm")) {
  #  surface = match.arg(surface)
  p  = cbind(obj$points$x, obj$points$y)
  p1 = p[obj$triangles[,1],]
  p2 = p[obj$triangles[,2],]
  p3 = p[obj$triangles[,3],]
  
  altitude = function(p1,p2,p3) {
    v = p1 - p2
    w = p3 - p2
    g = v - rowSums(v*w) / rowSums(w*w) * w
    g / rowSums(g*v)
  }
  v = p1 - p2
  w = p3 - p2
  a = abs(v[,2]*w[,1] - v[,1]*w[,2])/2
  area = sum(a)
  sapply(surface, function(surface) {
    h = obj$points[[surface]]
    h = cbind(h[obj$triangles[,1]], h[obj$triangles[,2]], h[obj$triangles[,3]])
    grad = h[,1] * altitude(p1,p2,p3) + h[,2] * altitude(p2,p3,p1) + h[,3] * altitude(p3,p1,p2)
    c(
      area = area,
      volume = sum(a*rowMeans(h)),
      min = min(h),
      max = max(h),
      mean = sum(a*rowMeans(h))/area,
      meanSquare = sum(a*rowMeans(h^2))/area, # todo: approximation
      meanGradSquare = sum(a * rowSums(grad^2))/area
    )
  })
}

#' Set the gap and offset for already created fracture
#' 
#' @param obj the fracture_geom object
#' @param gap the new gap
#' @param closed the probability of closed gap
#' @param offset the offset by which the fracture should be moved vertical
#' 
#' @export
set_gap = function(obj, gap, closed, offset) {
  if (missing(gap)) {
    if (missing(closed)) {
      gap = obj$gap
    } else {
      gap=-qnorm(closed, mean=0, sd=sqrt(obj$var.diff))
    }
  }
  if (missing(offset)) offset = obj$offset
  dgap = gap - obj$gap
  doffset = offset - obj$offset
  obj$points$f1 = obj$points$f1 + doffset + dgap/2
  obj$points$f2 = obj$points$f2 + doffset - dgap/2
  obj$points$fm = obj$points$fm + doffset
  obj$points$h = obj$points$h + dgap
  obj$gap = gap
  obj$offset = offset
  obj
}
