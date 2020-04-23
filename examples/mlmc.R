
tab$perm11

tab$level = factor(tab$refine)
levels(tab$level) = LETTERS[seq_len(nlevels(tab$level))]

reshape2::dcast(tab,seed~refine, value.var = "perm11")


library(Matrix)


MyCov = function(x,y) {
  x$i = seq_len(nrow(x))
  y$i = seq_len(nrow(y))
  xy = merge(x,y,by="seed")
  sparseMatrix(i=xy$i.x, j=xy$i.y, x = K[as.integer(xy$level.x) + 3*(as.integer(xy$level.y)-1)],dims=c(nrow(x),nrow(y)))
}
MyH = function(x) {
  sparseMatrix(i=seq_len(nrow(x)), j=as.integer(x$level), x = 1, dims=c(nrow(x),3))
}


K = diag(3)
Z = sparseMatrix(i=c(),j=c(),x=numeric(0),dims=c(3,3))

X = tab[,c("seed","level")]
X$val = tab$perm11/tab$perm_gap
Y = expand.grid(seed=unique(tab$seed), level=factor(levels(tab$level)))

for (it in seq_len(50)) {
M = MyCov(X,X)
N = MyH(X)
MAT = rbind(cbind(M,N),cbind(t(N),Z))
RET = solve(MAT, c(X$val,rep(0,3)))
nM = MyCov(Y,X)
nN = MyH(Y)
nMAT = cbind(nM,nN)
Y$val = as.vector(nMAT %*% RET)
Y$var = colSums(t(nMAT) * solve(MAT, t(nMAT)))
YW = reshape2::dcast(Y, seed~level, value.var = "val")
matplot(YW[,2], YW[,1:3+1],pch=16,cex=0.3)
r = sqrt(sum((K - cov(YW[,1:3+1]))^2))
print(r)
K = cov(YW[,1:3+1])
mu = RET[1:3+nrow(X)]
}

XW = reshape2::dcast(X, seed~level, value.var = "val")
matplot(XW[,3], XW[,1:3+1],pch=16,cex=0.3)
cov(XW[,1:3+1],use = "pairwise")
K
colMeans(XW[,1:3+1],na.rm = TRUE)
mu

rgl::plot3d(XW[,1:3+1])

nYW = merge(YW,data.frame(seed = XW$seed, lev = 3-rowSums(is.na(XW))),by="seed")

rgl::plot3d(nYW[,1:3+1],col=nYW$lev)

w = rep(0,nrow(X)+3)
w[nrow(X)+3] = 1

nx = data.frame(seed=0, level=1:3)

c = diag(MyCov(nx,nx))
b = rbind(MyCov(X,nx), t(MyH(nx)))
g = solve(MAT, b)
dim(g)

(w %*% g)^2 / (c - colSums(b*g))

w %*% solve(MAT) %*% w - c(w,0) %*% solve(cbind(rbind(MAT,t(b)),rbind(b,c))) %*% c(w,0)


w %*% solve(MAT) %*% w
eigen(solve(MAT)[nrow(X)+1:3,nrow(X)+1:3])

K[1]/nrow(X)
