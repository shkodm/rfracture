library(Matrix)

N = 5000
cost = c(1,2,4)
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

levs = factor(paste0("lev",1:4))
nlevs = nlevels(levs)
x = expand.grid(lev = levs,seed = 1:20)

fun = function(seed,lev) Ytable[seed+nrow(Ytable)*(as.integer(lev)-1)]


times = rep(0,nlevs)
evals = rep(0,nlevs)
#fun2 = function(seed,lev) Ytable[seed+nrow(Ytable)*(as.integer(lev)-1)]
fun2 = function(seed,lev) {
  #lev = c(2,1,3)[lev2]
  library(rfracture)
  alpha = 4.5
  sd = 0.001
  power.iso = function(f) sd^2*ifelse(f<5,1,(f/5)^-alpha)
  corr.profile = function(lambda) 0
  G = 12
  gap = 1/(2*G)
  refine = 2^(as.integer(lev))
  #refine = round(2^seq(2,4,len=4))[as.integer(lev)]
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

loglik = function(SK, SX) {
  eSK = eigen(SK)
  vSX = SX %*% eSK$vectors
  vSX = colSums(vSX^2)
  eval = eSK$values
  selval = eval > max(eval)*1e-6
  dval = ifelse(selval,1/eval,0)
  k = nrow(SX)
  -(k*sum(log(eval[selval]*2*pi)) + sum(dval*vSX))/2
}

TF = c(TRUE,FALSE)
TF = as.matrix(do.call(expand.grid,lapply(seq_len(nlevs),function(i) TF)))
TF = lapply(seq_len(nrow(TF)),function(i)as.vector(TF[i,]))

loglik_na = function(K,X) {
  sum(sapply(TF, function(sel) {
    inX = is.na(X)
    rsel = apply(cbind(!inX[,sel,drop=FALSE],inX[,!sel,drop=FALSE]),1,all)
    if (sum(rsel) > 0) {
      loglik(K[sel,sel,drop=FALSE], X[rsel,sel,drop=FALSE])
    } else {
      0
    }
  }))
}

p_to_K = function(p) {
  mt = matrix(0,nlevs,nlevs)
  mt[row(mt) >= col(mt)] = p
  mt %*% t(mt)
}


x$val = sapply(seq_len(nrow(x)),function(i) fun(x$seed[i],x$lev[i]))
x0 = x





Z = sparseZero(nlevs,nlevs)

K = diag(nrow=nlevs)
p = K[row(K) >= col(K)]
log = NULL
for (ad in seq_len(4000)) {
  #x = x[order(x$seed,x$lev),]
  
  #if (round(log2(ad),digits = 7) == log2(ad)) {
  if (nrow(x) %% 100 == 0) {
    x_tab = reshape2::dcast(x, seed~lev, value.var = "val")
    X = as.matrix(x_tab[,seq_len(nlevs)+1])
    X = apply(X,2,function(x) x - mean(x,na.rm=TRUE))
    ret = optim(p,function(p) -loglik_na(p_to_K(p),X),method="L-BFGS-B")
    p = ret$par
    K = p_to_K(p)
    print(K)
  }
#  K = realK
  
  xH = MyH(x)
  xS = MyCov(x,x)
  xM = rbind(cbind(xS,xH),cbind(t(xH),Z))
  MY = solve(xM,c(x$val,rep(0,nlevs)))
  mu = MY[nrow(x)+1:3]
  mlmc_mu = mean((x_tab[,2])[!is.na(x_tab[,2]) & is.na(x_tab[,3]) & is.na(x_tab[,4])]) +
    mean((x_tab[,3]-x_tab[,2])[!is.na(x_tab[,2]) & !is.na(x_tab[,3]) & is.na(x_tab[,4])]) +
    mean((x_tab[,4]-x_tab[,3])[!is.na(x_tab[,3]) & !is.na(x_tab[,4])])
  tmp_fun = function(x) var(x)/length(x)
  mlmc_var = tmp_fun((x_tab[,2])[!is.na(x_tab[,2]) & is.na(x_tab[,3]) & is.na(x_tab[,4])]) +
    tmp_fun((x_tab[,3]-x_tab[,2])[!is.na(x_tab[,2]) & !is.na(x_tab[,3]) & is.na(x_tab[,4])]) +
    tmp_fun((x_tab[,4]-x_tab[,3])[!is.na(x_tab[,3]) & !is.na(x_tab[,4])])
  x_tab = reshape2::dcast(x, seed~lev, value.var = "val")
  
h = rep(0,nrow(x)+nlevs)
h[nrow(x)+nlevs] = 1

x_tab = reshape2::dcast(x, seed~lev, value.var = "val")
x_na = is.na(x_tab[,1:nlevs+1])
x_na = rbind(x_na,TRUE)
rownames(x_na) = c(x_tab$seed,max(x_tab$seed)+1)

a = do.call(rbind,lapply(seq_len(nlevs), function(i) cbind(x_na,i)[x_na[,i],,drop=FALSE]))
a = unique(a)
nx = data.frame(seed=as.integer(rownames(a)),lev = levs[a[,nlevs+1]])

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
log = rbind(log,c(sum(cost[x$lev]),var,mlmc_var,mu,mlmc_mu,K[3,3],var(x$val[x$lev=="lev3"])))
#print(vardif)

sel = which.max(vardifpercost)
#print(a[sel,drop=FALSE])
#sel = sample(seq_along(vardifpercost),size = 1, prob = vardifpercost)

add_x = nx[sel,]
add_x$val = fun(add_x$seed,add_x$lev)
print(add_x)
x = rbind(x,add_x)
if (nrow(log) %% 10 == 0) {
  par(mfrow=c(2,2))
#  matplot(log[,1],log[,1:6+9],log="y",type="l",col=rep(1:3,times=2),lty=rep(1:2,each=3))
#  plot(-log[,2],log="xy")
#  abline(0,-1)
  #plot(log[,1],log[,17]); abline(h=2)
  #lines(log[,1],2 + sqrt(-log[,2]),lty=2)
  #lines(log[,1],2 - sqrt(-log[,2]),lty=2)
  fin = log[nrow(log),6]
  matplot(log[,1],log[,6:7],type="l",ylim=quantile(log[,6:7],probs = c(0.01,0.99),na.rm = TRUE),lty=1)
  abline(h=fin)
  lines(log[,1],fin + sqrt(log[,2]),lty=2)
  lines(log[,1],fin - sqrt(log[,2]),lty=2)
  lines(log[,1],fin + sqrt(log[,3]),lty=2,col=2)
  lines(log[,1],fin - sqrt(log[,3]),lty=2,col=2)
  
  matplot(log[,1],log[,8:9],type="l",ylim=quantile(log[,8:9],probs = c(0.1,0.9)),lty=1)
  
  ltop = K[nlevs,nlevs]/min(log[,1])*cost[nlevs]
  llow = min(log[,2:3],ltop,na.rm = TRUE)
  matplot(log[,1],log[,2:3],log="xy",type="l",lty=1,ylim=c(llow,ltop),xlim=c(min(log[,1]),cost[nlevs]*K[nlevs,nlevs]/llow))
  #lines(log2[,1],log2[,2],col=8)
  #lines(seq_len(nrow(log))*cost[3],realK[3,3]/seq_len(nrow(log)),lty=2,col=2)
  #lines(seq_len(nrow(log))*cost[3],K[3,3]/seq_len(nrow(log)),lty=2,col=1)
  for (i in seq_len(nlevs)) abline(log10(K[i,i]*cost[i]),-1,lty=ifelse(i == nlevs,2,3))
  #abline(log10(K[3,3]*cost[3]),-1,lty=2)
  #abline(log10(K[2,2]*cost[2]),-1,lty=3)
  #abline(log10(K[1,1]*cost[1]),-1,lty=3)
  wy = table(x$lev)
  wx = barplot(wy)
  text(wx,wy,wy, adj=c(0.5,-0.3) )
}
}



#plot(-log[,2],log="xy")
#abline(0,-1)


