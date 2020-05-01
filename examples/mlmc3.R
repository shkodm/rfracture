library(Matrix)
library(future)

plan(multisession)

levs = factor(paste0("lev",1:4))
target = 4
N = 5000
cost = c(1,2,4)
#cost = c(1,1,1)
set.seed(123)
seeds = sample.int(.Machine$integer.max,size = N, replace = TRUE)

nlevs = nlevels(levs)

fun = function(seed,lev) Ytable[seed+nrow(Ytable)*(as.integer(lev)-1)]


times = rep(0,nlevs)
evals = rep(0,nlevs)
#fun2 = function(seed,lev) Ytable[seed+nrow(Ytable)*(as.integer(lev)-1)]
fun = function(seed,lev) {
  fun2 = function(seed,lev) {
    #lev = c(2,1,3)[lev2]
    library(rfracture)
    alpha = 4.5
    sd = 0.001
    power.iso = function(f) sd^2*ifelse(f<5,1,(f/5)^-alpha)
    corr.profile = function(lambda) 0
    G = 12
    gap = 1/(2*G)
    refine = 2^(as.integer(lev)+1)
    #refine = round(2^seq(2,4,len=4))[as.integer(lev)]
    ret = fracture_geom(refine=refine, corr.profile = corr.profile, power.iso = power.iso, seed=seed)
    ret2 = set_gap(ret, gap = gap)
    ret2 = slice(ret2,  value="above")
    ret2 = cut(ret2)
    res = solve_reynolds(ret2, method = "direct")
    res$perm[1,1]/res$perm_gap
  }
  lev = as.integer(lev)
  tm = system.time({ret = fun2(seed,lev); gc();})
  list(time=as.numeric(tm[1]), value=ret)
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
TF = TF[sapply(TF, any)]

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


x = expand.grid(lev = levs,seed = 1:6)
x$val = NA

fut = NULL

Z = sparseZero(nlevs,nlevs)

K = diag(nrow=nlevs)
p = K[row(K) >= col(K)]
log = NULL
plot_time = Sys.time()
K_time = Sys.time()
K_n = 0
for (ad in seq_len(4000)) {
  #x = x[order(x$seed,x$lev),]

  to_calc = setdiff(which(is.na(x$val)),sapply(fut, function(x) x$index))
  if (length(to_calc) >= 1) {
    cat("Adding",length(to_calc),"futures to run...\n")
    fut = c(
      lapply(to_calc,function(i) list(index=i, value=future(fun(seeds[x$seed[i]],x$lev[i])))),
      fut
    )
  }
  
  #if (round(log2(ad),digits = 7) == log2(ad)) {
#  if (difftime(Sys.time(), K_time, units = "sec") > 5) {
  sel = sapply(fut, function(x) resolved(x$value))
  if (length(sel) == 0) sel=logical()
  if (K_n == 0) sel = rep(TRUE, length(fut))
  if (sum(sel) > 0) {
    cat("Evaluating",sum(sel),"futures...\n")
    for (f in fut[sel]) {
      idx = f$index
      lev = x$lev[idx]
      val = value(f$value)
      x$val[idx] = val$value
      times[lev] = times[lev] + val$time
      evals[lev] = evals[lev] + 1
    }
    fut = fut[!sel]
    cost = times/evals
  }
  
  if (sum(!is.na(x$val)) - K_n >= 20) {
    x_tab = reshape2::dcast(x, seed~lev, value.var = "val")
    X = as.matrix(x_tab[,seq_len(nlevs)+1])
    X = apply(X,2,function(x) x - mean(x,na.rm=TRUE))
    ret = optim(p,function(p) -loglik_na(p_to_K(p),X),method="L-BFGS-B")
    p = ret$par
    K = p_to_K(p)
    print(K)
    K_time = Sys.time()
    K_n = sum(!is.na(x$val))
  }
#  K = realK
  
  xH = MyH(x)
  xS = MyCov(x,x)
  xM = rbind(cbind(xS,xH),cbind(t(xH),Z))
  MY = solve(xM,c(x$val,rep(0,nlevs)))
  mu = MY[nrow(x)+seq_len(nlevs)]

  h = rep(0,nrow(x)+target)
  h[nrow(x)+nlevs] = 1

  x_tab = reshape2::dcast(x, seed~lev, value.var = "seed")
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

  var = -as.vector(h %*% solve(xM,h))
  vardif = as.vector((h %*% g)^2 / (c - colSums(g*b)))
  vardifpercost = vardif / cost[nx$lev]
  
  log = rbind(log,c(
    sum(cost[x$lev]), # cost
    var,              # var of result
    mu[target],       # mu
    K[target,target], # variance
    var(x$val[x$lev==levs[target]]), # variance from fine grid
    tapply(vardifpercost,nx$lev,max)
  ))

  #sel = which.max(vardifpercost)
  #print(a[sel,drop=FALSE])
  sel = sample(seq_along(vardifpercost),size = 1, prob = vardifpercost)
  
  add_x = nx[sel,]
  #add_x$val = fun(add_x$seed,add_x$lev)$value
  add_x$val = NA
  print(add_x)
  x = rbind(x,add_x)
  if (difftime(Sys.time(), plot_time, units = "sec") > 3) {
    plot_time = Sys.time()
    par(mfrow=c(2,2))
    fin = log[nrow(log),3]
    plot(log[,1],log[,3],type="p",ylim=quantile(log[,3],probs = c(0.01,0.99),na.rm = TRUE),lty=1)
    abline(h=fin)
    lines(log[,1],fin + sqrt(log[,2]),lty=2)
    lines(log[,1],fin - sqrt(log[,2]),lty=2)
    lines(log[,1],fin + sqrt(log[,3]),lty=2,col=2)
    lines(log[,1],fin - sqrt(log[,3]),lty=2,col=2)
    
    #matplot(log[,1],log[,4:5,drop=FALSE],type="p",ylim=quantile(log[,4:5],probs = c(0.1,0.9),na.rm = TRUE),lty=1)
    #matplot(log[,1:nlevs+5],log="y",type="l",lty=1)
    matplot(log[,1:nlevs+5]/rowSums(log[,1:nlevs+5,drop=FALSE]),ylim=c(0,1),type="l",lty=1)
    ltop = K[target,target]/min(log[,1])*cost[target]
    llow = min(log[,2],ltop,na.rm = TRUE)
    plot(log[,1],log[,2],log="xy",type="l",lty=1,ylim=c(llow,ltop),xlim=c(min(log[,1]),cost[nlevs]*K[nlevs,nlevs]/llow))
    for (i in seq_len(nlevs)) abline(log10(K[i,i]*cost[i]),-1,lty=ifelse(i == target,2,3))
    wy = table(x$lev)
    wx = barplot(wy)
    text(wx,wy,wy, adj=c(0.5,-0.3) )
  }
}



#plot(-log[,2],log="xy")
#abline(0,-1)


