library(rfracture)
alpha = 4.5
sd = 0.001
#sd = 0.000
power.iso = function(f) sd^2*ifelse(f<5,1,(f/5)^-alpha)
corr.profile = function(lambda) ifelse(lambda<0.5,0.1,0.9)
corr.profile = function(lambda) 0.2

G = 12
mar = 1/G/32 * 0.6

refine = 32
gap_n = 12

ret = fracture_geom(refine=refine, corr.profile = corr.profile, power.iso = power.iso, seed=123)
  
gap = gap_n / (G*32)
ret2 = set_gap(ret, gap = gap)
ret2 = slice(ret2,  value="above")
ret2 = slice(ret2,  level=1/G - 2*mar, flatten="above")
if (require(rgl)) fracture3d(cut(ret2),edge.col = 4, vertex.col = 4)
obj = cut(ret2)
    
#' Summary of surfaces
#' 
#' @param obj fracture_geom object
#' @param surface the variables over which to gather statistics
#' 
#' @export

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
  g1 = altitude(p1,p2,p3)
  g2 = altitude(p2,p3,p1)
  g3 = altitude(p3,p1,p2)
  A11 = rowSums(g1*g1)*a
  A12 = rowSums(g1*g2)*a
  A13 = rowSums(g1*g3)*a
  A22 = rowSums(g2*g2)*a
  A23 = rowSums(g2*g3)*a
  A33 = rowSums(g3*g3)*a
  t1 = obj$triangles[,1]
  t2 = obj$triangles[,2]
  t3 = obj$triangles[,3]
  M = data.frame(i = c( t1, t1, t1, t2, t2, t2, t3, t3, t3),
                 j = c( t1, t2, t3, t1, t2, t3, t1, t2, t3),
                 v = c(A11,A12,A13,A12,A22,A23,A13,A23,A33))
  library(Matrix)
  M = rbind(M[M$i > 1, ], data.frame(i=1, j=1:nrow(obj$points), v=1))
  M = sparseMatrix(i=M$i, j=M$j, x=M$v)
  
  plot(as.vector(M %*% rep(1,nrow(M))))
  sum(as.vector(M %*% obj$points$y) * obj$points$y)
  sum(a)
  
  ret = solve(M, c(0,rep(1, nrow(M)-1)))

  plot(ret)
  