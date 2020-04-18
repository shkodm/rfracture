
estimate_wishart = function(res, dof, gam, ...) {
  if (missing(dof)) dof = 10
  if (missing(gam)) gam = 1
  lGamma_2 = function(n) log(pi)/2 + lgamma(n) + lgamma(n-1/2)
  tr = function(x) sum(diag(x))
  sld = sum(sapply(res, function(x) determinant(x,logarithm = TRUE)$modulus))
  str = sum(sapply(res, function(x) tr(x)))
  N = length(res)
  loglik = function(n,gam) -N*n*p/2*log(2*gam) - N*lGamma_2(n/2) + (n-p-1)/2*sld -str/(2*gam)
  ret = optim(c(dof,gam), function(par) -loglik(par[1],par[2]), ...)
  ret$dof = ret$par[1]
  ret$gam = ret$par[2]
  ret$loglik = -ret$value
  ret
}

estimate_wishart_lam = function(res, dof, gam, ...) {
  if (missing(dof)) dof = 10
  if (missing(gam)) gam = 1
  lGamma_2 = function(n) log(pi)/2 + lgamma(n) + lgamma(n-1/2)
  lam = t(sapply(res, function(x) eigen(x)$value))
  sld = sum(log(lam))
  str = sum(lam)
  N = length(res)
  loglik = function(n,gam) -N*n*p/2*log(2*gam) - N*lGamma_2(n/2) + (n-p-1)/2*sld - str/(2*gam)
  ret = optim(c(dof,gam), function(par) -loglik(par[1],par[2]), ...)
  ret$dof = ret$par[1]
  ret$gam = ret$par[2]
  ret$loglik = -ret$value
  ret
}

estimate_wishart_lam2 = function(res, dof, ...) {
  if (missing(dof)) dof = 10
  lGamma_2 = function(n) log(pi)/2 + lgamma(n) + lgamma(n-1/2)
  lam = t(sapply(res, function(x) eigen(x)$value))
  sld = sum(log(lam))
  str = sum(lam)
  N = length(res)
  loglik = function(n) {
    g = N*n*p/2
    -g*(log(str) - log(g) + 1) - N*lGamma_2(n/2) + (n-p-1)/2*sld
  }
  ret = optim(c(dof), function(par) -loglik(par[1]), ...)
  ret$dof = ret$par[1]
  ret$gam = str/(N*p*ret$dof)
  ret$loglik = -ret$value
  ret
}


estimate_wishart_lam_mean = function(res, dof, mean, ...) {
  if (missing(dof)) dof = 10
  if (missing(mean)) mean = 0
  lGamma_2 = function(n) log(pi)/2 + lgamma(n) + lgamma(n-1/2)
  lam = t(sapply(res, function(x) eigen(x)$value))
  N = length(res)
  loglik = function(n,mean) {
    g = N*n*p/2
    -g*(log(sum(lam-mean)) - log(g) + 1) - N*lGamma_2(n/2) + (n-p-1)/2 * sum(log(lam-mean))
  }
  ret = optim(c(dof,mean), function(par) -loglik(par[1],par[2]), ...)
  ret$dof = ret$par[1]
  ret$mean = ret$par[2]
  ret$gam = sum(lam-ret$mean)/(N*p*ret$dof)
  ret$loglik = -ret$value
  ret
}


ret1 = estimate_wishart(res)
ret2 = estimate_wishart_lam(res)
ret3 = estimate_wishart_lam2(res)
ret4 = estimate_wishart_lam_mean(res, dof=20, mean=0)

lam = t(sapply(res, function(x) eigen(x)$value))
plot(lam)
abline(v=ret3$gam*ret3$dof)

M = sapply(res, function(x) c(x[1,1],x[1,2],x[2,1],x[2,2]))
matrix(rowMeans(M),2,2)

ret3$gam*ret3$dof
ret4$gam*ret4$dof + ret4$mean

estimate.wishart_m = function(res, dof, gam, mean, ...) {
  if (missing(dof)) dof = 10
  if (missing(gam)) gam = 1
  if (missing(mean)) mean = 0
  lGamma_2 = function(n) log(pi)/2 + lgamma(n) + lgamma(n-1/2)
  tr = function(x) sum(diag(x))
  sld = sum(sapply(res, function(x) determinant(x,logarithm = TRUE)$modulus))
  str = sum(sapply(res, function(x) tr(x)))
  N = length(res)
  loglik = function(n,gam) -N*n*p/2*log(2*gam) - N*lGamma_2(n/2) + (n-p-1)/2*sld -str/(2*gam)
  ret = optim(c(dof,gam), function(par) -loglik(par[1],par[2]), ...)
  ret$dof = ret$par[1]
  ret$gam = ret$par[2]
  ret$loglik = -ret$value
  ret
}

