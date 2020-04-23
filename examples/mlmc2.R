library(Matrix)

N = 300
set.seed(124)
Ytable = matrix(rnorm(N*3),N,3,byrow = TRUE)
Ytable = cbind(Ytable,1) %*% matrix(c(1/3,1/10,1,1,0,1/10,1,1.5,0,0,1,2),4,3)

#rgl::plot3d(Ytable)

levs = factor(paste0("lev",1:3))
nlevs = nlevels(levs)
x = expand.grid(seed = 1:4, lev = levs)

fun = function(seed,lev) Ytable[seed+nrow(Ytable)*(as.integer(lev)-1)]

x$val = fun(x$seed,x$lev)

x = x[seq_len(nrow(x)-1),]

x_tab = reshape2::dcast(x, seed~lev, value.var = "val")

sparseZero = function(n,m) sparseMatrix(i=integer(0),j=integer(0),x=numeric(0),dims=c(n,m))

MyCov = function(x,y) {
  X = data.frame(seed=x$seed, lev=as.integer(x$lev), i=seq_len(nrow(x)))
  Y = data.frame(seed=y$seed, lev=as.integer(y$lev), i=seq_len(nrow(y)))
  
  XY = merge(X,Y,by="seed")
  sparseMatrix(
    i=XY$i.x,
    j=XY$i.y,
    x=K[XY$lev.x + (XY$lev.y-1)*nlevs],
    dims=c(nrow(X),nrow(Y)))
}

MyH = function(x) {
  sparseMatrix(
    i=seq_len(nrow(x)),
    j=as.integer(x$lev),
    x=1,
    dims=c(nrow(x),nlevs))
}

Z = sparseZero(3,3)

K = diag(nrow=nlevs)

for (it in seq_len(50)) {
  xH = MyH(x)
  xS = MyCov(x,x)
  xM = rbind(cbind(xS,xH),cbind(t(xH),Z))
  
  xe = expand.grid(seed = unique(x$seed), lev = levs)
  xe$val = cbind(MyCov(xe,x),MyH(xe)) %*% solve(xM,c(x$val,rep(0,nlevs)))
  
  xe_tab = reshape2::dcast(xe, seed~lev, value.var = "val")
  K_ = cov(xe_tab[,seq_len(nlevs)+1])
  r = sqrt(sum((K-K_)^2))
  print(r)
  if (r < 1e-10) break;
  K = K_
}

h = rep(0,nrow(x)+3)
h[nrow(x)+3] = 1

nx = expand.grid(seed = max(x$seed)+1, lev = levs)
nx = rbind(nx,data.frame(seed=4,lev="lev3"))
nS = MyCov(x,nx)
nH = MyH(nx)
b = rbind(nS,t(nH))
g = solve(xM,b)
c = diag(MyCov(nx,nx))

vardif = (h %*% g)^2 / (c - colSums(g*b))

cost = c(1,2,3)

vardif / cost[nx$lev]
x

sum(h*solve(xM,h))
sum(c(h,0)*solve(cbind(rbind(xM,b[,3]),c(b[,3],c[3])),c(h,0)))


