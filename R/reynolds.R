#' Calculated the premability tensor by solving a Reynolds lubrication equation in the fracture
#' @param obj the fracture_geom object
#' @param q.fun the dependence of the flux on the heights, by default h^3/12
#' @importFrom Matrix sparseMatrix solve t 
#' @export
solve_reynolds = function(obj, q.fun=function(h) h^3/12) {
  
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
  
  i = rep(0, nrow(obj$points))
  i[sel] = as.integer(factor(paste(round(obj$points$x[sel],5) %% 1,round(obj$points$y[sel],5) %% 1,sep="_")))
  i[!sel] = seq_len(sum(!sel)) + max(i)
  
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
  t1 = i[obj$triangles[,1]]
  t2 = i[obj$triangles[,2]]
  t3 = i[obj$triangles[,3]]
  
  M = data.frame(i = c( t1, t1, t1, t2, t2, t2, t3, t3, t3),
                 j = c( t1, t2, t3, t1, t2, t3, t1, t2, t3),
                 v = c(A11,A12,A13,A12,A22,A23,A13,A23,A33))
  
  A = sparseMatrix(i=c(t1,t2,t3), j=rep(1,length(t1)*3), x=c(a,a,a)/3)
  M = sparseMatrix(i=M$i, j=M$j, x=M$v)
  
  BM = rbind(cbind(M,A),cbind(t(A),1))
  
  D = sparseMatrix(i=c(t1,t2,t3,t1,t2,t3),j=rep(1:2,each=length(t1)*3), x=c(g1[,1],g2[,1],g3[,1],g1[,2],g2[,2],g3[,2])*trinta)
  
  RHS = rbind(D,0)  
  
  #The hard part
  X = solve(BM, RHS)
  
  VOL = sum((h1+h2+h3)/3*a)
  perm_diff = -as.matrix(t(X) %*% RHS)
  perm_hom = sum(trinta)
  I = diag(nrow=2)
  list(
    perm_diff = perm_diff,
    perm_hom = I * perm_hom,
    perm = perm_diff + I * perm_hom,
    avgu = (perm_diff + I * perm_hom) / VOL,
    perm_gap = I * q.fun(obj$gap),
    avgu_gap = I * q.fun(obj$gap) / obj$gap,
    avgu_h = I * perm_hom / VOL
  ) 
}
