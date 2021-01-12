#' Calculated the premability tensor by solving a Reynolds lubrication equation in the fracture
#' @param obj the fracture_geom object
#' @param q.fun the dependence of the flux on the heights, by default h^3/12
#' @importFrom Matrix sparseMatrix solve t 
#' @export
solve_reynolds = function(obj, q.fun=function(h) h^3/12, method=c("direct","iterative","system")) {
  method = match.arg(method)
  p  = cbind(obj$points$x, obj$points$y)
  p1 = p[obj$triangles[,1],]
  p2 = p[obj$triangles[,2],]
  p3 = p[obj$triangles[,3],]
  
  h1 = obj$points$h[obj$triangles[,1]]
  h2 = obj$points$h[obj$triangles[,2]]
  h3 = obj$points$h[obj$triangles[,3]]
  
  trint = (q.fun(h1)+q.fun(h2)+q.fun(h3))/3
  
  altitude = function(p1,p2,p3) {
    v = p1 - p2
    w = p3 - p2
    g = v - rowSums(v*w) / rowSums(w*w) * w
    g / rowSums(g*v)
  }

  sel = FALSE
  sel = sel | abs(obj$points$x - 0) < 1e-12
  sel = sel | abs(obj$points$x - obj$width) < 1e-12
  sel = sel | abs(obj$points$y - 0) < 1e-12
  sel = sel | abs(obj$points$y - obj$width) < 1e-12

  v = p1 - p2
  w = p3 - p2
  a = abs(v[,2]*w[,1] - v[,1]*w[,2])/2
  trinta = trint * a
  g1 = altitude(p1,p2,p3)
  g2 = altitude(p2,p3,p1)
  g3 = altitude(p3,p1,p2)
  A11 = rowSums(g1*g1)*trinta
  A12 = rowSums(g1*g2)*trinta
  A13 = rowSums(g1*g3)*trinta
  A22 = rowSums(g2*g2)*trinta
  A23 = rowSums(g2*g3)*trinta
  A33 = rowSums(g3*g3)*trinta
  
  i = rep(0, nrow(obj$points))
  i[sel] = as.integer(factor(paste(round(obj$points$x[sel],5) %% 1,round(obj$points$y[sel],5) %% 1,sep="_")))
  i[!sel] = seq_len(sum(!sel)) + max(i)
  
  sel = trinta >= max(trinta) * 1e-6
  
  t1 = i[obj$triangles[,1]]
  t2 = i[obj$triangles[,2]]
  t3 = i[obj$triangles[,3]]
  
  M = data.frame(i = c( t1, t1, t1, t2, t2, t2, t3, t3, t3),
                 j = c( t1, t2, t3, t1, t2, t3, t1, t2, t3),
                 v = c(A11,A12,A13,A12,A22,A23,A13,A23,A33))
  
  A = sparseMatrix(i=c(t1,t2,t3), j=rep(1,length(t1)*3), x=c(a,a,a)/3)
  M = sparseMatrix(i=M$i, j=M$j, x=M$v)
  
  BM = rbind(cbind(M,A),cbind(t(A),0))
  
  D = sparseMatrix(i=c(t1,t2,t3,t1,t2,t3),j=rep(1:2,each=length(t1)*3), x=c(g1[,1],g2[,1],g3[,1],g1[,2],g2[,2],g3[,2])*trinta)
  
  RHS = rbind(D,0)  

  #The hard part
  if (method == "direct") {
    X = try(cbind(
      solve(BM, RHS[,1]),
      solve(BM, RHS[,2])
    ))
    if (inherits(X, "try-error")) {
      errors = "Solve failed"
      perm_diff = matrix(NA,2,2)
    } else {
      perm_diff = -as.matrix(t(X) %*% RHS)
      errors = NA
    }
    iter = NA
  } else if (method == "iterative") {
    ret = Rlinsolve::lsolve.gmres(BM,RHS)
    X = ret$x
    iter = ret$iter
    errors = ret$errors
    perm_diff = -as.matrix(t(X) %*% RHS)
  } else if (method == "system") {
    return (list(matrix= BM, rhs=RHS))
  }
  
  VOL = sum((h1+h2+h3)/3*a)
  area = sum(a)
  hmean = VOL/area
  perm_hom = sum(trinta)
  I = diag(nrow=2)
  list(
    perm_diff = perm_diff,
    perm_hom = perm_hom,
    perm = perm_diff + I * perm_hom,
    perm_gap = q.fun(obj$gap)*(obj$width*obj$width),
    volume = VOL,
    area = area,
    hmean = hmean,
    perm_hmean = q.fun(hmean)*area,
    iter = iter,
    errors = errors,
    lambdas = as.vector(X[nrow(D)+1,])
  ) 
}

#' Simple Conjugate Gradient implementation
#' 
#' @param A Matrix
#' @param b Right-hand-side
#' @param itmax Maximum number of iterations
#' @param epsjump The level of residual drop needed
#' @param eps Absolute level of residual needed, by default equal to epsjump times initial residual
#' @param x0 Initial guess on solution
#' 
#' @export
solve_cg = function(A, b, itmax = 800, epsjump=1e-6, eps, x0 = rep(0, ncol(A))) {
  x = x0
  r = b - A %*% x
  q = r
  r2 = sum(r^2)
  res = r2
  if (missing(eps)) eps = epsjump * r2
  for (i in 1:itmax) {
    r2 = sum(r^2)
    res = c(res, r2)
    if (r2 < eps) break
    Aq = A %*% q
    alpha = sum(r * q)/sum(q * Aq)
    x = x + alpha*q
    r = r - alpha*Aq
    beta = sum(r * Aq)/sum(q * Aq)
    q = r - beta * q
  }
  list(x=x, iterations=i, residual=r2, eps=eps, res=res)
}

