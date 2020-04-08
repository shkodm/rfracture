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
    
  p  = cbind(obj$points$x, obj$points$y)
  p1 = p[obj$triangles[,1],]
  p2 = p[obj$triangles[,2],]
  p3 = p[obj$triangles[,3],]

  h1 = obj$points$h[obj$triangles[,1]]
  h2 = obj$points$h[obj$triangles[,2]]
  h3 = obj$points$h[obj$triangles[,3]]
  
  trint = (fun(h1)+fun(h2)+fun(h3))/3
  trinta = trint * a
  
  altitude = function(p1,p2,p3) {
    v = p1 - p2
    w = p3 - p2
    g = v - rowSums(v*w) / rowSums(w*w) * w
    g / rowSums(g*v)
  }
  
  
  sel = FALSE
  sel = sel | abs(obj$points$x - 0) < 1e-12
  sel = sel | abs(obj$points$x - 1) < 1e-12
  sel = sel | abs(obj$points$y - 0) < 1e-12
  sel = sel | abs(obj$points$y - 1) < 1e-12
  
  i = rep(0, nrow(obj$points))
  i[sel] = as.integer(factor(paste(round(obj$points$x[sel],5) %% 1,round(obj$points$y[sel],5) %% 1,sep="_")))
  i[!sel] = seq_len(sum(!sel)) + max(i)
  
  
  v = p1 - p2
  w = p3 - p2
  a = abs(v[,2]*w[,1] - v[,1]*w[,2])/2
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
  library(Matrix)
  A = sparseMatrix(i=c(t1,t2,t3), j=rep(1,length(t1)*3), x=c(a,a,a)/3)
  M = sparseMatrix(i=M$i, j=M$j, x=M$v)
  
  BM = rbind(cbind(M,A),cbind(t(A),1))
  
  D = sparseMatrix(i=c(t1,t2,t3,t1,t2,t3),j=rep(1:2,each=length(t1)*3), x=c(g1[,1],g2[,1],g3[,1],g1[,2],g2[,2],g3[,2])*trinta)

  RHS = rbind(D,0)  
  dim(RHS)
  dim(BM)
  X = solve(BM, RHS)
  dim(X)
  plot3d(obj$points$x, obj$points$y, X[i,1] - obj$points$x)
  plot3d(obj$points$x, obj$points$y, X[i,2] - obj$points$y)
  
  dim(X)
  dim(RHS)
  VOL = sum((h1+h2+h3)/3*a)
  RES = t(X) %*% RHS / VOL
  
  RES
  eigen(RES)
  
  plot(as.vector(M %*% rep(1,nrow(M))))
  sum(as.vector(M %*% obj$points$y) * obj$points$y)
  sum(a)
  
  ret = solve(M, c(0,rep(1, nrow(M)-1)))

  plot(ret)
  

  
  
  i = seq_len(nrow(obj$points))
  sel = abs(obj$points$x - 0) < 1e-12
  a1 = data.frame(i = i[sel],v = round(obj$points$y[sel],5))
  sel = abs(obj$points$x - 1) < 1e-12
  a2 = data.frame(j = i[sel],v = round(obj$points$y[sel],5))
  if (nrow(a1) != nrow(a2)) stop("wrong lengths")
  a3 = merge(a1,a2)
  if (nrow(a3) != nrow(a1)) stop("wrong length of merge")
  
  i[a3$i] = i[a3$j]

sel = abs(obj$points$y - 0) < 1e-12
a1 = data.frame(i = i[sel],v = round(obj$points$x[sel],5))
sel = abs(obj$points$y - 1) < 1e-12
a2 = data.frame(j = i[sel],v = round(obj$points$x[sel],5))
if (nrow(a1) != nrow(a2)) stop("wrong lengths")
a3 = merge(a1,a2)
if (nrow(a3) != nrow(a1)) stop("wrong length of merge")

i[a3$i] = i[a3$j]

plot(i)

plot3d(obj$points$x, obj$points$y, i)





plot(i)
range(i)


h1 = obj$points$h[obj$triangles[,1]]
h2 = obj$points$h[obj$triangles[,2]]
h3 = obj$points$h[obj$triangles[,3]]

fun = function(h) h^3/12
plot(fun((h1+h2+h3)/3)/((fun(h1)+fun(h2)+fun(h3))/3))


surface_summary(obj)

tr = as.vector(t(obj$triangles))
plot3d(obj$points$x, obj$points$y, (X[i,1] - obj$points$x))
triangles3d(obj$points$x[tr], obj$points$y[tr], (X[i,1] - obj$points$x)[tr],col=3)

