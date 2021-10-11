#' Calculated the premability tensor by solving a Reynolds lubrication equation in the fracture
#' @param obj the fracture_geom object
#' @param q.fun the dependence of the flux on the heights, by default h^3/12
#' @importFrom Matrix sparseMatrix solve t 
#' @export
solve_reynolds = function(obj, q.fun=function(h) h^3/12, method=c("direct","iterative","system"), rule, alt=0, hmult='none', ...) {
  method = match.arg(method)
  if (is.character(alt)) {
    if (alt %in% names(obj$points)) {
      alt = obj$points$fm
    } else {
      stop("Unknown 'alt' value in solve_reynolds, possible options:", paste(names(obj$points,sep=",")))
    }
  }
  if (is.numeric(alt)) {
    if (! length(alt) %in% c(1,nrow(obj$points))) {
      stop("Wrong length of 'alt' in solve_reynolds.")
    }
  } else {
    stop("'alt' in solve_reynolds should be either field name or vector of values")
  }
  p  = cbind(obj$points$x, obj$points$y, alt)
  p1 = p[obj$triangles[,1],]
  p2 = p[obj$triangles[,2],]
  p3 = p[obj$triangles[,3],]
  
  h1 = obj$points$h[obj$triangles[,1]]
  h2 = obj$points$h[obj$triangles[,2]]
  h3 = obj$points$h[obj$triangles[,3]]
  
  grad = function(p1,p2,p3) {
    V2 = p1 - p2
    V1 = p3 - p2
    M11 = rowSums(V1*V1)
    M12 = rowSums(V1*V2)
    M22 = rowSums(V2*V2)
    DET = M11*M22 - M12*M12
    (-M12*V1 + M11*V2) / DET
  }
  
  sel = FALSE
  sel = sel | abs(obj$points$x - 0) < 1e-12
  sel = sel | abs(obj$points$x - obj$width) < 1e-12
  sel = sel | abs(obj$points$y - 0) < 1e-12
  sel = sel | abs(obj$points$y - obj$width) < 1e-12
  
  v = p1 - p2
  w = p3 - p2
  
  n = cbind(v[,2]*w[,3] - v[,3]*w[,2],v[,3]*w[,1] - v[,1]*w[,3],v[,1]*w[,2] - v[,2]*w[,1])
  a = sqrt(rowSums(n*n))
  n = n / a
  a = a / 2
  
  if (is.character(hmult)) {
    if (hmult == "slant") {
      hmult = abs(n[,3])
    } else if (hmult == "none") {
      hmult = 1
    }
  }
  if (is.numeric(hmult)) {
    if (! length(hmult) %in% c(1, nrow(obj$triangles))) {
      stop("Wrong length of 'hmult' in solve_reynolds")
    }
  } else {
    stop("'hmult' in solve_reynolds, should be 'slant', 'none' or vector of values")
  }
  trint = (q.fun(hmult*h1)+q.fun(hmult*h2)+q.fun(hmult*h3))/3
  
  trinta = trint * a
  g1 = grad(p1,p2,p3)
  g2 = grad(p2,p3,p1)
  g3 = grad(p3,p1,p2)
  vMw = function(v,w,m) m*rowSums(v*w)
  A11 = vMw(g1,g1,trinta)
  A12 = vMw(g1,g2,trinta)
  A13 = vMw(g1,g3,trinta)
  A22 = vMw(g2,g2,trinta)
  A23 = vMw(g2,g3,trinta)
  A33 = vMw(g3,g3,trinta)
  
  i = rep(0, nrow(obj$points))
  i[sel] = as.integer(factor(paste(round(obj$points$x[sel]/obj$width,5) %% 1,round(obj$points$y[sel]/obj$width,5) %% 1,sep="_")))
  i[!sel] = seq_len(sum(!sel)) + max(i)
  
  t1 = obj$triangles[,1]
  t2 = obj$triangles[,2]
  t3 = obj$triangles[,3]
  
  Per = sparseMatrix(i=i, j=seq_along(i), x=1)
  M = sparseMatrix( i = c( t1, t1, t1, t2, t2, t2, t3, t3, t3),
                    j = c( t1, t2, t3, t1, t2, t3, t1, t2, t3),
                    x = c(A11,A12,A13,A12,A22,A23,A13,A23,A33) )
  A = sparseMatrix(i=c(t1,t2,t3), j=rep(1,length(t1)*3), x=c(a,a,a)/3)
  D = M %*% p[,1:2]
  
  BM = rbind(cbind(Per %*% M %*% t(Per),Per %*% A),cbind(t(Per %*% A),0))
  RHS = rbind(Per %*% D,0)  

  errors = NULL
  # The hard part
  if (method == "direct") {
    X = try(cbind(
      solve(BM, RHS[,1]),
      solve(BM, RHS[,2])
    ))
    if (inherits(X, "try-error")) {
      errors = "Solve failed"
      X = rep(NA,prod(dim(RHS)))
      dim(X) = dim(RHS)
    }
    iter = NA
  } else if (method == "iterative") {
    ret = solve_cg(BM, RHS, ...)
    X = ret$x
    iter = ret$iterations
    errors = ret$errors
  } else if (method == "system") {
    return (list(matrix= BM, rhs=RHS))
  }
  
  perm_diff = - as.matrix(t(RHS) %*% X)
  perm_hom = as.matrix(t(p[,1:2]) %*% M %*% p[,1:2])
  
  VOL = sum((h1+h2+h3)/3*a)
  area = sum(a)
  base_area = obj$width^2
  hmean = VOL/area
  perm_diff = perm_diff / base_area
  perm_hom  = perm_hom  / base_area
  list(
    perm_diff = perm_diff,
    perm_hom = perm_hom,
    perm = perm_diff + perm_hom,
    perm_gap = q.fun(obj$gap),
    perm_hmean = q.fun(hmean),
    volume = VOL,
    area = area,
    base_area = base_area,
    hmean = hmean,
    iter = iter,
    errors = errors,
    lambdas = as.vector(X[nrow(X),])
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
solve_cg = function(A, b, itmax = 800, epsjump=1e-6, eps, x0 = rep(0, ncol(A)), print.every=500) {
  x = x0
  r = b - A %*% x
  q = r
  r2 = sum(r^2)
  res = r2
  if (missing(eps)) eps = epsjump * r2
  for (i in 1:itmax) {
    r2 = sum(r^2)
    res = c(res, r2)
    if ((i %% print.every) == 0) cat("CG: iteration:",i," res:", r2, "\n",sep="")
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

