#' Generates a sequence of numbers between 0 and 1, without 1
#' @export
seq_1 = function(n) (seq_len(n)-1)/n

#' Generates a circular sequence of frequencies (0 to n/2 and -n/2 to 0)
#' @export
seq_circ = function(n) { x = 1:n-1; ifelse(x > n/2, x-n, x) }

#' Expands a grid of sequences (outputs a matrix)
#' @export
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
ordered_rnorm_spectrum = function(f, k=2, seed, length_one=FALSE, fracture_coeffs_generator=NULL) {
  if (! missing(seed)) set.seed(seed)
  f = apply(f,2,round)
  if (inherits(f,"matrix")) f = as.data.frame(f)
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
  if(is.null(fracture_coeffs_generator)) {
        # Original random generation
        RN = matrix(rnorm(totsize * k * 2), 2, k*totsize)
    } else {
        # Get values from generator for specified fracture
        RN = fracture_coeffs_generator(seed, totsize, k)
    }
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
#' power.spectrum = function(f) f^{-2.5}
#' par(mfrow=c(2,2))
#' ret = fracture_matrix(dims=c(50,50),span = diag(2), power.iso = power.spectrum, seed=seed)
#' plot(ret)
#' ret = fracture_matrix(dims=c(50,50),span = matrix(c(1,1,-1,1),2,2), power.iso = power.spectrum, seed=seed)
#' plot(ret)
#' ret = fracture_matrix(dims=c(50,50),span = diag(2)*2, power.iso = power.spectrum, seed=seed)
#' plot(ret)
#' ret = fracture_matrix(dims=c(50,50),span = diag(2)*2,period=diag(2)*2, power.iso = power.spectrum, seed=seed)
#' plot(ret)
#' 
#' @import stats
#' @export
fracture_matrix = function(
    dims = c(10,10),
    span = diag(nrow=length(dims)),
    period = diag(nrow=length(dims)),
    shift = rep(0,length(dims)),
    power.iso,
    power.spectrum = function(f) power.iso(sqrt(rowSums(f*f))),
    cov.iso,
    cov.structure = function(d) cov.iso(sqrt(rowSums(d*d))),
    corr.profile = function(k) 0,
    closed = 0.1, gap, seed, corr.method = c("midline","top","mixed"), length_one = FALSE, gauss.order = 1, cov.offsets.number = 5,
    fracture_coeffs_generator = NULL
    ) {
  span = as.matrix(span)
  period = as.matrix(period)
  corr.method = match.arg(corr.method)
  n = length(dims)
  p_ = expand_seq(dims,seq_1)
  p = p_ %*% span
  f_ = expand_seq(dims,seq_circ)
  f = f_ %*% solve(t(span))


  power = rep(0,nrow(f))
  
  
  f_per = f %*% t(period)
  sel = rowSums(abs(round(f_per) - f_per) > 1e-6) == 0
  coef = matrix(0, nrow(f_per), 2)
  coef[sel,] = ordered_rnorm_spectrum(f_per[sel,,drop=FALSE], k = 2, seed = seed, length_one = length_one, fracture_coeffs_generator = fracture_coeffs_generator)

  do.cov   = ! (missing(cov.iso) && missing(cov.structure))
  do.power = ! (missing(power.iso) && missing(power.spectrum))
  if (do.cov) {
    if (do.power) stop("supplied both coveriance structure and power spectrum")
    offsets = expand_seq(rep(cov.offsets.number*2,length(dims)),seq_circ)
    cov.field = 0
    for (i in 1:nrow(offsets)) {
      d = p - as.matrix(offsets[rep(i,nrow(p)),]) %*% period
      c = cov.structure(d)
      cov.field = cov.field + c
    }
#    d = ((p %*% solve(period) + 0.5) %% 1 - 0.5) %*% period
#    cov.field = cov.structure(d)
    dim(cov.field) = dims
#    image(cov.field)
    
    power = fft(cov.field)
    power = power / prod(dims)
    print(c(range(Re(power)),range(Im(power))))
    power = Re(power)
    power = ifelse(power > 0, power, 0)
    power = as.vector(power)
    power[!sel] = 0
    power_mult = 1
    power = power * power_mult
  } else if (do.power) {
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
    
    power_mult = 1/det(period)
    power = power * power_mult
  } else {
    stop("neither covariance structure nor power spectrum supplied")
  }  
  freq = sqrt(rowSums(f^2))
  wavelength = 1/freq
  
  #Sn = n*pi^(n/2)/gamma(n/2+1)/2
  
  rad = sqrt(power)
  corr.prof = corr.profile(wavelength)
  fshift = f %*% shift * pi * 2
  if (any(abs(corr.prof) > 1)) stop("Correlation outside of [-1,1] interval")
  if (corr.method == "midline") { # midline always the same
    # ang = acos(corr)/2
    # M = matrix(c(cos(ang),sin(ang),cos(ang),-sin(ang)),2,2)
    corr.angle = acos(corr.prof)/2
    M11 =   cos(corr.angle) * rad * exp(- 0.5i * fshift);
    M12 =   sin(corr.angle) * rad * exp(- 0.5i * fshift);
    M21 =   cos(corr.angle) * rad * exp(  0.5i * fshift);
    M22 = - sin(corr.angle) * rad * exp(  0.5i * fshift);
  } else if (corr.method == "mixed") { # nice mix of two random variables
    # ang = asin(corr)/2
    # M = matrix(c(cos(ang),sin(ang),sin(ang),cos(ang)),2,2)
    corr.angle = asin(corr.prof)/2
    corr.coef = list((cos(corr.angle) * coef[,1] + sin(corr.angle) * coef[,2]) * rad * exp(- 0.5i * fshift),
                     (sin(corr.angle) * coef[,1] + cos(corr.angle) * coef[,2]) * rad * exp(  0.5i * fshift))
    M11 =   cos(corr.angle) * rad * exp(- 0.5i * fshift);
    M12 =   sin(corr.angle) * rad * exp(- 0.5i * fshift);
    M21 =   sin(corr.angle) * rad * exp(  0.5i * fshift);
    M22 =   cos(corr.angle) * rad * exp(  0.5i * fshift);
  } else if (corr.method == "top") { # top surface always the same
    # ang = asin(corr)
    # M = matrix(c(1,0,sin(ang),cos(ang)),2,2)
    corr.angle = asin(corr.prof)
    corr.coef = list((coef[,1]) * rad,
                     (sin(corr.angle) * coef[,1] + cos(corr.angle) * coef[,2]) * rad * exp(  1i * fshift))
    M11 =                 1 * rad * 1;
    M12 =                 0 * rad * 1;
    M21 =   sin(corr.angle) * rad * exp(    1i * fshift);
    M22 =   cos(corr.angle) * rad * exp(    1i * fshift);
  } else {
    stop("unknown corr.method")
  }
  
  corr.coef = list(M11 * coef[,1] + M12 * coef[,2],
                   M21 * coef[,1] + M22 * coef[,2])

  fields = lapply(corr.coef, function(x) {
    dim(x) = dims
    Re(fft(x, inverse = TRUE))
  })
  
  ret = list()
  c11 = Re(sum(M11 * Conj(M11) + M12 * Conj(M12)))
  c12 = Re(sum(M11 * Conj(M21) + M12 * Conj(M22)))
  c21 = Re(sum(M21 * Conj(M11) + M22 * Conj(M12)))
  c22 = Re(sum(M21 * Conj(M21) + M22 * Conj(M22)))
  cov.theoretical = matrix(c(c11,c12,c21,c22),2,2)
  cov.final = cov(cbind(as.vector(fields[[1]]),as.vector(fields[[2]])))
  # E(((f1+f2)/2)^2) = E(f1^2 + f2^2 + 2*f1*f2)/4
  var.midline = Re((c11+c12+c21+22)/4)
  # E((f1-f2)^2) = E(f1^2 + f2^2 - 2*f1*f2)
  var.diff = Re(c11-c12-c21+c22)
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
    offset = 0,
    prob.closed = pnorm(-gap,0,sd=sqrt(var.diff)),
    length_one = length_one,
    gauss.order = gauss.order,
    spec1D = tapply(power,abs(f[,1]),sum)/2,
    power_mult = power_mult
  )
  class(ret) = "fracture_matrix"
  ret
}


#' Calculate Matern covariance function and its power spectrum
#' 
#' @param r distance
#' @param nu smoothness parameter
#' @param l scale
#' @export
matern = function(r, nu, l=1) {
  a = sqrt(2*nu)*r/l
  ifelse(a > 0,
    2^(1-nu)/gamma(nu)*(a)^nu*besselK(a,nu = nu),
    1
  )
}

#' @rdname matern
#' @export
matern.f = function(s, nu, l=1, D) {
  (2*sqrt(pi))^D*gamma(nu+D/2)*(2*nu)^nu / (gamma(nu)*l^(2*nu))*(2*nu/l^2 + 4*pi^2*s^2)^(-(nu+D/2))
}


