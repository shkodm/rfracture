library(Matrix)

N = 5000
cost = c(1,2,3)
#cost = c(1,1,1)
set.seed(123)

seeds = sample.int(.Machine$integer.max,size = N, replace = TRUE)
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
times = rep(0,nlevs)
evals = rep(0,nlevs)
#fun2 = function(seed,lev) Ytable[seed+nrow(Ytable)*(as.integer(lev)-1)]
fun2 = function(seed,lev) {
  alpha = 4.5
  sd = 0.001
  power.iso = function(f) sd^2*ifelse(f<5,1,(f/5)^-alpha)
  corr.profile = function(lambda) 0
  G = 12
  gap = 1/(2*G)
  refine = 2^(as.integer(lev)+1)
  ret = fracture_geom(refine=refine, corr.profile = corr.profile, power.iso = power.iso, seed=seed)
  ret2 = set_gap(ret, gap = gap)
  ret2 = slice(ret2,  value="above")
  ret2 = cut(ret2)
  res = solve_reynolds(ret2, method = "direct")
  res$perm[1,1]/res$perm_gap
}
fun = function(seed,lev) {
  lev = as.integer(lev)
  tm = system.time({ret = fun2(seeds[seed],lev)})
  times[lev] <<- times[lev] + tm[1]
  evals[lev] <<- evals[lev] + 1
  cost <<- times/evals
  ret
}

x$val = sapply(seq_len(nrow(x)),function(i) fun(x$seed[i],x$lev[i]))

#x = x[seq_len(nrow(x)-1),]

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
  #K = realK
  x_tab = reshape2::dcast(x, seed~lev, value.var = "val")
  # cov_biased = function(m) {
  #   m = as.matrix(m)
  #   m = apply(m,2,function(x) x-mean(x,na.rm=TRUE))
  #   nv = colSums(!is.na(m))
  #   m1 = m; m1[ is.na(m1)] = 0;
  #   m2 = m; m2[!is.na(m2)] = 1; m2[is.na(m2)] = 0;
  #   N2 = sqrt(outer(nv,nv))
  #   N1 = t(m2) %*% m2
  #   (t(m1) %*% m1)/N2
  # }
  # K = cov_biased(x_tab[,1:3+1])
  vars = c(
    var(x_tab[,2],na.rm=TRUE),
    var(x_tab[,3]-x_tab[,2],na.rm=TRUE),
    var(x_tab[,4]-x_tab[,3],na.rm=TRUE)
  )
  mt = matrix(c(1,0,0,-1,1,0,0,-1,1),3,3)
  #K = solve(t(mt)) %*% diag(vars) %*% solve(mt)
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
  #break;
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

#sel = which.max(vardifpercost)
sel = sample(1:6, size = 1, prob = vardifpercost)

add_x = nx[sel,]
add_x$val = fun(add_x$seed,add_x$lev)
print(add_x)
x = rbind(x,add_x)
if (nrow(log) %% 10 == 0) {
  par(mfrow=c(2,2))
  matplot(log[,1],log[,1:6+9],log="y",type="l",col=rep(1:3,times=2),lty=rep(1:2,each=3))
#  plot(-log[,2],log="xy")
#  abline(0,-1)
  fin = log[nrow(log),17]
  plot(log[,1],log[,17],type="l"); abline(h=fin)
  lines(log[,1],fin + sqrt(log[,2]),lty=2)
  lines(log[,1],fin - sqrt(log[,2]),lty=2)
  mean(x$val[x$lev=="lev3"])
  plot(log[,1],log[,2],log="xy",type="l",lty=1)
  lines(seq_len(nrow(log))*cost[3],K[3,3]/seq_len(nrow(log)),lty=2)
  barplot(table(x$lev))
}
}


#plot(-log[,2],log="xy")
#abline(0,-1)


#x_tab = reshape2::dcast(x, seed~lev, value.var = "val")

