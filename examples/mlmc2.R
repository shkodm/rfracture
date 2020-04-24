library(Matrix)

N = 5000
cost = c(1,2,3)
#cost = c(1,1,1)
set.seed(123)
Ytable = matrix(rnorm(N*3),N,3,byrow = TRUE)
Ymat = matrix(c(1/3,1/13,1,1,0,1/10,1,1.5,0,0,1,2),4,3)
#Ymat = matrix(c(1,0,0,1,0,1,0,1.5,0,0,1,2),4,3)
Ytable = cbind(Ytable,1) %*% Ymat
Ytable = cbind(Ytable,1) %*% matrix(c(1,0,0,1,0,1,0,1.5,0,0,1,2),4,3)

realK = t(Ymat) %*% diag(c(1,1,1,0)) %*% Ymat

#realK = diag(3)
#rgl::plot3d(Ytable)

levs = factor(paste0("lev",1:3))
nlevs = nlevels(levs)
x = expand.grid(lev = levs,seed = 1:20)

fun = function(seed,lev) Ytable[seed+nrow(Ytable)*(as.integer(lev)-1)]

x$val = fun(x$seed,x$lev)

x = x[seq_len(nrow(x)-1),]

x_tab = reshape2::dcast(x, seed~lev, value.var = "val")

sparseZero = function(n,m) sparseMatrix(i=integer(0),j=integer(0),x=numeric(0),dims=c(n,m))

MyCov = function(x,y,BASE=K) {
  X = data.frame(seed=x$seed, lev=as.integer(x$lev), i=seq_len(nrow(x)))
  Y = data.frame(seed=y$seed, lev=as.integer(y$lev), i=seq_len(nrow(y)))
  
  XY = merge(X,Y,by="seed")
  sparseMatrix(
    i=XY$i.x,
    j=XY$i.y,
    x=BASE[XY$lev.x + (XY$lev.y-1)*nlevs],
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
log = NULL
for (ad in seq_len(2000)) {
  x = x[order(x$seed,x$lev),]
  K = realK
  #x_tab = reshape2::dcast(x, seed~lev, value.var = "val")
  #K = cov(x_tab[,1:3+1],use="pairwise") 
for (it in seq_len(50)) {
  xH = MyH(x)
  xS = MyCov(x,x)
  xM = rbind(cbind(xS,xH),cbind(t(xH),Z))
  MY = solve(xM,c(x$val,rep(0,nlevs)))
  mu = MY[nrow(x)+1:3]
  xe = expand.grid(seed = unique(x$seed), lev = levs)
  xe$val = cbind(MyCov(xe,x),MyH(xe)) %*% MY
  
  xe_tab = reshape2::dcast(xe, seed~lev, value.var = "val")
  K_ = cov(xe_tab[,seq_len(nlevs)+1])
  r = sqrt(sum((K-K_)^2))
  #print(r)
  break;
  if (r < 1e-10) break;
  K = K_
}

h = rep(0,nrow(x)+3)
h[nrow(x)+3] = 1

nx = expand.grid(seed = max(x$seed)+1, lev = levs)

lev_seed = sapply(levs, function(l) max(x$seed[x$lev == l])) + 1
nx = rbind(
  expand.grid(seed = max(x$seed) + 1, lev=levs),
  data.frame(seed = sapply(levs, function(l) max(x$seed[x$lev == l]))+1, lev = levs)
)
#nx = unique(nx)

nS = MyCov(x,nx)
nH = MyH(nx)
b = rbind(nS,t(nH))
g = solve(xM,b)
c = diag(MyCov(nx,nx))

mh = solve(xM,h)
var = -as.vector(h %*% solve(xM,h))
wei = as.vector(mh)[seq_len(nrow(x))]
length(wei)
var2 = as.vector(wei %*% MyCov(x,x,realK) %*% wei)
vardif = as.vector((h %*% g)^2 / (c - colSums(g*b)))

vardifpercost = vardif / cost[nx$lev]
log = rbind(log,c(sum(cost[x$lev]),var,var2,vardif,vardifpercost, mu))
#print(vardif)

sel = which.max(vardifpercost)

add_x = nx[sel,]
add_x$val = fun(add_x$seed,add_x$lev)
print(add_x)
x = rbind(x,add_x)
if (nrow(log) %% 10 == 0) {
  par(mfrow=c(2,2))
  matplot(log[,1],log[,1:6+8],log="y",type="l",col=rep(1:3,times=2),lty=rep(1:2,each=3))
#  plot(-log[,2],log="xy")
#  abline(0,-1)
  #plot(log[,1],log[,17]); abline(h=2)
  #lines(log[,1],2 + sqrt(-log[,2]),lty=2)
  #lines(log[,1],2 - sqrt(-log[,2]),lty=2)
  matplot(log[,2:3],log="y",type="l",lty=1)
  lines(realK[3,3]/seq_len(nrow(log)),lty=2)
  barplot(table(x$lev))
}
}


#plot(-log[,2],log="xy")
#abline(0,-1)


#x_tab = reshape2::dcast(x, seed~lev, value.var = "val")

