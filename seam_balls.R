library(fields)

par(mfrow=c(1,1))

P = obj$points
i = obj$triangles
iv = as.vector(t(i))

source("~/Dropbox/rtri.R")
T = data.frame(i1=i[,1],i2=i[,2],i3=i[,3])
h = P$h[i]
dim(h) = dim(i)
T$h = rowMeans(h)

sel =        P$x[T$i1] > 0 & P$x[T$i2] > 0 & P$x[T$i3] > 0
sel = sel & (P$y[T$i1] > 0 & P$y[T$i2] > 0 & P$y[T$i3] > 0)
sel = sel & (P$x[T$i1] < 1 & P$x[T$i2] < 1 & P$x[T$i3] < 1)
sel = sel & (P$y[T$i1] < 1 & P$y[T$i2] < 1 & P$y[T$i3] < 1)
T = T[sel,]

set.seed(100)
K = 4000
Rmax = 0.01
Rmin = 0.01
margin = 0.003
margin_opt = margin/3
B = data.frame()
ndel = 100
dlog = NULL
for (iterations in 1:600) {
  if (K+ndel > nrow(B)) {
    #k = min((K-nrow(B)),200)
    k = min(K+ndel - nrow(B),200)
    p = pmax(0, T$h - 2 * (Rmin+margin))
    
    nB = data.frame(tri=sample.int(nrow(T), k, prob = p,replace = TRUE), r=runif(k,Rmin,Rmax))
    nB$h1 = P$h[T$i1[nB$tri]] - 2 * (nB$r+margin)
    nB$h2 = P$h[T$i2[nB$tri]] - 2 * (nB$r+margin)
    nB$h3 = P$h[T$i3[nB$tri]] - 2 * (nB$r+margin)
    sel = nB$h1+nB$h2+nB$h3 > 0
    nB = nB[sel,,drop=FALSE]
    ret = rtri(nrow(nB),nB$h1,nB$h2,nB$h3)
    nB$y = P$x[T$i1[nB$tri]]*ret$w1 + P$x[T$i2[nB$tri]]*ret$w2 + P$x[T$i3[nB$tri]]*ret$w3
    nB$z = P$y[T$i1[nB$tri]]*ret$w1 + P$y[T$i2[nB$tri]]*ret$w2 + P$y[T$i3[nB$tri]]*ret$w3
    nB$f1 = P$f1[T$i1[nB$tri]]*ret$w1 + P$f1[T$i2[nB$tri]]*ret$w2 + P$f1[T$i3[nB$tri]]*ret$w3
    nB$f2 = P$f2[T$i1[nB$tri]]*ret$w1 + P$f2[T$i2[nB$tri]]*ret$w2 + P$f2[T$i3[nB$tri]]*ret$w3
    nB$x = nB$f2 + nB$r + margin + ret$h
    if (! (all(nB$x > nB$f2 + nB$r + margin) && all(nB$x < nB$f1 - nB$r - margin))) stop("balls do not fit in triangles")
    B = rbind(B, nB)
    cat("Created", nrow(nB), "balls\n")
  }

  kdel = max(10, nrow(B) - K)

  X = as.matrix(B[, c("x", "y", "z")])
  clear3d()
  points3d(X,col=1)
  for (i in 2:10) {
    X2 = rbind(X,as.matrix(P[,c("f1","x","y")]),as.matrix(P[,c("f2","x","y")]))
    ds = fields.rdist.near(X2,X,delta = 2*(Rmax+margin_opt))
    tds = data.frame(v = ds$ra, i = ds$ind[,1], j = ds$ind[,2])
    tds$v = tds$v - ifelse(tds$i > nrow(X),0,B$r[tds$i]) - B$r[tds$j]
    sel = tds$i != tds$j
    tds = tds[sel,]
    print(range(tds$v))
    #  plot(sort(tds$v))
#    par(mfrow=c(2,1))
#    plot(tds$i,tds$j,col=ifelse(tds$v<0,2,0),cex=0.1,xlim=c(1,nrow(X)),ylim=c(1,nrow(X)))
#    plot(tds$i,tds$j,col=ifelse(tds$v<0,2,0),cex=0.1,xlim=c(1,nrow(X2)),ylim=c(1,nrow(X)))
    #  o = order(tds$v)[1:(kdel*4)]
    #  tds = tds[o,]
    #  sel = tds$v < 0
    #  tds = tds[sel,]
    p = X2[tds$i,] - X[tds$j,]
    pl = sqrt(rowSums(p^2))
    p = p / pl * pmax(0,margin - tds$v)
    p = sapply(1:3, function(k) tapply(p[,k],tds$j,mean))
    pi = tapply(tds$j,tds$j,function(x) x[1])
    #clear3d()
#    mysegments3d(X[pi,],X[pi,]-p,col=i)
    X[pi,] = X[pi,] - p
#    points3d(X,col=i)
  }
  
  B[, c("x", "y", "z")] = X
  ds = fields.rdist.near(X,delta = 2*(Rmax+margin_opt))
  tds = data.frame(v = ds$ra, i = ds$ind[,1], j = ds$ind[,2])
  tds$v = tds$v - B$r[tds$i] - B$r[tds$j]
  sel = tds$i != tds$j
  tds = tds[sel,]
  o = order(tds$v)
  tds = tds[o,]
  tds = tds[tds$v < 0,]
  delsel = NULL
  while(nrow(tds) > 0) {
    sel = tds$i[1]
    delsel = c(delsel, sel)
    tds = tds[tds$i != sel & tds$j != sel, , drop = FALSE]
  }
  cat("Deleted", length(delsel), "balls\n")
  dlog = rbind(dlog,data.frame(n=nrow(B)))
  plot(dlog$n)
  if (length(delsel) > 0)  B = B[-delsel,]
  if (nrow(B) >= K) {
    cat("Finished.\n")
    break;
  }
}

iv = as.vector(t(obj$triangles))
col = colorRamp(c("black", "red", "yellow"))((h - min(h)) / (max(h) -
                                                               min(h)))
col = rgb(col, max = 255)
clear3d()
plot3d(diag(3))
clear3d()
triangles3d(P$f1[iv], P$x[iv], P$y[iv], col = rep(col, each = 3))
spheres3d(B$x, B$y, B$z, col = "green", radius = B$r)
#triangles3d(P$f2[iv], P$x[iv], P$y[iv], col = 1,alpha=0.7)


names(B)
write.csv(B[,c("x","y","z","r")],file="seam_balls.csv",row.names=FALSE)


# atom-ID atom-type diameter density x y z
atoms = data.frame(id=1:nrow(B), type=1, diam=B$r*2, dens=2, x=B$x+0.1, y=B$y, z=B$z)
f = file("seam.gran_shift.data",open = "w")
cat("LIGGGHTS Description\n",file=f)
cat("\n",file=f)
cat(nrow(B),"atoms\n",file=f)
cat("1 atom types\n",file=f)
cat("0 bond types\n",file=f)
cat("0 angle types\n",file=f)
cat("\n",file=f)
cat(" 0.0 0.2 xlo xhi\n",file=f)
cat(" 0.0 1.0 ylo yhi\n",file=f)
cat(" 0.0 1.0 zlo zhi\n",file=f)
cat("\n",file=f)
cat("Atoms\n",file=f)
cat("\n",file=f)
write.table(atoms,file=f,row.names=FALSE,col.names=FALSE)
close(f)

range(atoms$x)


# 
# 
# X = as.matrix(B[, c("x", "y", "z")])
# 
# ds = fields.rdist.near(X,delta = 2*(Rmax+margin_opt))
# tds = data.frame(v = ds$ra, i = ds$ind[,1], j = ds$ind[,2])
# tds$v = tds$v - B$r[tds$i] - B$r[tds$j]
# sel = tds$i != tds$j
# tds = tds[sel,]
# 
# p = X[tds$i,] - X[tds$j,]
# pl = sqrt(rowSums(p^2))
# r2 = B$r[tds$i] + B$r[tds$j]
# p = p / pl * pmax(r2-pl,0)
# p = sapply(1:3, function(k) tapply(p[,k],tds$i,sum))
# pi = tapply(tds$i,tds$i,function(x) x[1])
# dim(p)
# dim(X[pi,])
# X[pi,] = X[pi,] - p
# 
# p = X[tds$i,] - X[tds$j,]
# pl2 = sqrt(rowSums(p^2))
# matplot((cbind(r2,pl2,pl)/r2)[order(pl/r2),],cex=0.2,pch=16)
# 
# B[, c("x", "y", "z")] = X
# 
# col = colorRamp(c("black", "red", "yellow"))((h - min(h)) / (max(h) -
#                                                                min(h)))
# col = rgb(col, max = 255)
# clear3d()
# plot3d(diag(3))
# clear3d()
# triangles3d(P$f1[iv], P$x[iv], P$y[iv], col = rep(col, each = 3))
# spheres3d(B$x, B$y, B$z, col = "green", radius = B$r)


B = B2[1:300,]

clear3d()
X = as.matrix(B[, c("x", "y", "z")])
points3d(X,col=1)
for (i in 2:10) {
  X2 = rbind(X,as.matrix(P[,c("f1","x","y")]),as.matrix(P[,c("f2","x","y")]))
  ds = fields.rdist.near(X2,X,delta = 2*(Rmax+margin_opt))
  tds = data.frame(v = ds$ra, i = ds$ind[,1], j = ds$ind[,2])
  tds$v = tds$v - ifelse(tds$i > nrow(X),0,B$r[tds$i]) - B$r[tds$j]
  sel = tds$i != tds$j
  tds = tds[sel,]
  print(range(tds$v))
#  plot(sort(tds$v))
  par(mfrow=c(2,1))
  plot(tds$i,tds$j,col=ifelse(tds$v<0,2,0),cex=0.1,xlim=c(1,nrow(X)),ylim=c(1,nrow(X)))
  plot(tds$i,tds$j,col=ifelse(tds$v<0,2,0),cex=0.1,xlim=c(1,nrow(X2)),ylim=c(1,nrow(X)))
#  o = order(tds$v)[1:(kdel*4)]
#  tds = tds[o,]
#  sel = tds$v < 0
#  tds = tds[sel,]
  p = X2[tds$i,] - X[tds$j,]
  pl = sqrt(rowSums(p^2))
  p = p / pl * pmax(0,margin - tds$v)
  p = sapply(1:3, function(k) tapply(p[,k],tds$j,mean))
  pi = tapply(tds$j,tds$j,function(x) x[1])
  #clear3d()
  mysegments3d(X[pi,],X[pi,]-p,col=i)
  X[pi,] = X[pi,] - p
  points3d(X,col=i)
}

points3d(p)
B[, c("x", "y", "z")] = X

clear3d()
triangles3d(P$f1[iv], P$x[iv], P$y[iv], col = rep(col, each = 3))
spheres3d(B$x, B$y, B$z, col = "green", radius = B$r)


clear3d()
points3d(X[pi,])
points3d(X[pi,]+p,col=3)
mysegments3d = function(a,b,...) {
  x = rbind(a,b)
  n = nrow(a)
  i = as.matrix(t(matrix(1:(2*n),n,2)))
  segments3d(x[i,],...)
}

clear3d()
points3d(X[pi,])
mysegments3d(X[pi,],X[pi,]-p,col=3)

dim(p)
dim(X[pi,])
