#' Estimate parameters of a isotropic Wishart distribution
#' 
#' @param res list of 2x2 matrices
#' @param dof initial guess for degrees of freedom
#' @param ... other parameters passed to optim
#' 
#' For Y with Wishart distribution with dof degrees of freedom and gam*I covariance function
#' Estimates dof and gam from sample.
#' 
#' @export
estimate_wishart = function(res, dof, ...) {
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
  ret$mean = 0
  ret
}

#' Estimate parameters of a isotropic Wishart distribution with shift
#' 
#' @rdname estimate_wishart
#' 
#' @export
estimate_wishart_mean = function(res, dof, mean, ...) {
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
