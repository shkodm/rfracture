# Generates a sequence of numbers between 0 and 1, without 1
seq_1 = function(n) (seq_len(n)-1)/n

# Generates a circular sequence of frequencies (0 to n/2 and -n/2 to 0)
seq_circ = function(n) { x = 1:n-1; ifelse(x > n/2, x-n, x) }

# Expands a grid of sequences (outputs a matrix)
expand_seq = function(x, fun=seq_len, matrix=TRUE) {
  x = do.call(expand.grid,lapply(x,fun))
  if (matrix) as.matrix(x) else x
}

#' Generate a matrix of random complex numbers in a consistent order
#' 
#' @param f table of frequencies (wave numbers)
#' @param k number of independent matrices
#' @param seed random seed
#' @param length_one if TRUE, return complex numbers of length 1
#' @examples
#' f = expand.grid(0:1,0:1)
#' ordered_rnorm_spectrum(f, seed=123)
#' @import stats
#' @export
ordered_rnorm_spectrum = function(f, k=2, seed, length_one=FALSE) {
  if (! missing(seed)) set.seed(seed)
  f = apply(f,2,round)
  if (class(f) == "matrix") f = as.data.frame(f)
  size = max(sapply(f,max))
  D = ncol(f)
  ftot = expand_seq(rep(size,D),function(n) -n:n,matrix=FALSE)
  nm = paste0("F",seq_len(D))
  names(f) = nm
  names(ftot) = nm
  fmax = do.call(pmax, lapply(ftot,abs))
  ord = do.call(order,c(list(fmax), ftot))
  ftot = ftot[ord,,drop=FALSE]
  totsize = nrow(ftot)
  
  ftot_op = -ftot
  f$i = seq_len(nrow(f))
  ftot$j = seq_len(nrow(ftot))
  ftot_op$k = seq_len(nrow(ftot_op))
  fmerge = merge(merge(ftot_op,ftot),f)
  RN = matrix(rnorm(totsize * k * 2), 2, k*totsize)
  RN = RN[1,] + 1i * RN[2,]
  RN = t(matrix(RN, k, totsize))
  
  if (length_one) {
    W = RN[fmerge$j,,drop=FALSE] + Conj(RN[fmerge$k,,drop=FALSE])
    W = W / Mod(W)
    #W[fmerge$j == 1,,drop=FALSE] = 0
  } else {
    W = RN[fmerge$j,,drop=FALSE]
  }
  ret = matrix(0,nrow(f),k)
  ret[fmerge$i,] = W
  ret
}

#' Generate a fracture field matrix
#' 
#' @param dims Dimensions of the fracture matrix (length N)
#' @param span base spanning the fracture (matrix NxN)
#' @param period base of periodicity parallelogram
#' @param power.spectrum power spectrum of the fields (function of frequency)
#' @param corr.profile correlation profile between (function of wave length)
#' @param closed the probability of fields touching
#' @param gap mean gap between upper and lower fields (overwrites closed)
#' @param seed random seed
#' @param length_one use complex numbers of length one for field generation
#' @references
#' Brown, S. R. (1995). Simple mathematical model of a rough fracture. Journal of Geophysical Research: Solid Earth, 100(B4), 5941-5952.
#' @examples
#' seed = 123
#' par(mfrow=c(2,2))
#' ret = fracture_matrix(dims=c(50,50),span = diag(2),seed=seed)
#' plot(ret)
#' ret = fracture_matrix(dims=c(50,50),span = matrix(c(1,1,-1,1),2,2),seed=seed)
#' plot(ret)
#' ret = fracture_matrix(dims=c(50,50),span = diag(2)*2,seed=seed)
#' plot(ret)
#' ret = fracture_matrix(dims=c(50,50),span = diag(2)*2,period=diag(2)*2,seed=seed)
#' plot(ret)
#' 
#' @import stats
#' @export
fracture_matrix = function(
    dims = c(10,10),
    span = diag(nrow=length(dims)),
    period = diag(nrow=length(dims)),
    power.iso,
    power.spectrum = function(f) power.iso(sqrt(rowSums(f*f))),
    corr.profile = function(k) 0,
    closed = 0.1, gap, seed, length_one = FALSE, gauss.order = 1) {
  n = length(dims)
  p_ = expand_seq(dims,seq_1)
  p = p_ %*% span
  f_ = expand_seq(dims,seq_circ)
  f = f_ %*% solve(t(span))


  power = rep(0,nrow(f))
  
  
  f_per = f %*% t(period)
  sel = rowSums(abs(round(f_per) - f_per) > 1e-6) == 0
  coef = matrix(0, nrow(f_per), 2)
  coef[sel,] = ordered_rnorm_spectrum(f_per[sel,,drop=FALSE], k = 2, seed = seed, length_one = length_one)

  if (gauss.order == 1) {
    power[sel] = power.spectrum(f_per[sel,] %*% solve(t(period)))
  } else {
    q = statmod::gauss.quad(gauss.order)
    qx = as.matrix(do.call(expand.grid,rep(list(q$nodes/2),n)))
    qw = do.call(expand.grid,rep(list(q$weights/2),n))
    qw = apply(qw,1,prod)
    selpow = 0
    self = t(f_per[sel,])
    for (i in seq_len(nrow(qx))) {
      selpow = selpow + qw[i]*power.spectrum(t(self + qx[i,]) %*% solve(t(period)))
    }
    power[sel] = selpow
  }
  power[1] = 0
  if (any(power < 0)) stop("Negative power spectrum")
  
  
  freq = sqrt(rowSums(f^2))
  wavelength = 1/freq
  
  #Sn = n*pi^(n/2)/gamma(n/2+1)/2
  
  rad = sqrt(power)
  corr.prof = corr.profile(wavelength)
  if (any(abs(corr.prof) > 1)) stop("Correlation outside of [-1,1] interval")
  corr.angle = atan(corr.prof)

  corr.coef = list((cos(corr.angle) * coef[,1] + sin(corr.angle) * coef[,2]) * rad,
    (sin(corr.angle) * coef[,1] + cos(corr.angle) * coef[,2]) * rad)

  fields = lapply(corr.coef, function(x) {
    dim(x) = dims
    Re(fft(x, inverse = TRUE))
  })
  
  ret = list()
  c1 = sum((power)[sel])
  c2 = sum((power * (2*sin(corr.angle)*cos(corr.angle)))[sel])
  cov.theoretical = matrix(c(c1,c2,c2,c1),2,2)
  cov.final = cov(cbind(as.vector(fields[[1]]),as.vector(fields[[2]])))
  var.midline = sum(((cos(corr.angle) + sin(corr.angle))^2*power)[sel])/2
  var.diff = sum(((cos(corr.angle) - sin(corr.angle))^2*power)[sel])*2
  var.prime = sum((power*(2*pi*freq)^2)[sel])
  if (missing(gap)) {
    if (!missing(closed)) {
      gap = -qnorm(closed,mean=0,sd=sqrt(var.diff))
    } else{
      gap = 0
    }
  }
  ret = list(
    points = p,
    f1 = fields[[1]] + gap/2,
    f2 = fields[[2]] - gap/2,
    dims = dims,
    span = span,
    period = period,
    cov.theoretical = cov.theoretical,
    cov.final = cov.final,
    var.midline = var.midline,
    var.diff = var.diff,
    var.prime = var.prime,
    power.spectrum = power.spectrum,
    corr.profile = corr.profile,
    gap = gap,
    prob.closed = pnorm(-gap,0,sd=sqrt(var.diff)),
    length_one = length_one
  )
  class(ret) = "fracture_matrix"
  ret
}
