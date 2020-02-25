#'Generates a sequence of numbers between 0 and 1, without 1
#'
#' @param n length of the sequence
seq_1 = function(n) (seq_len(n)-1)/n

#'Generates a circular sequence of frequencies (0 to n/2 and -n/2 to 0)
#'
#' @param n length of the sequence
seq_circ = function(n) { x = 1:n-1; ifelse(x > n/2, x-n, x) }

#'Expands a grid of sequences (outputs a matrix)
#'
#' @param x lengths of sequences to use
expand_seq = function(x, fun=seq_len, matrix=TRUE) {
  x = do.call(expand.grid,lapply(x,fun))
  if (matrix) as.matrix(x) else x
}

#' Generate a matrix of random complex numbers in a consistent order
#' 
#' @param f
#' @param k number of independent matrices
#' @param seed random seed
#' @param length_one if TRUE, return complex numbers of length 1
#' @example 
#' f = expand.grid(0:1,0:1)
#' ordered_rnorm_spectrum(f, seed=123)
#' @export
ordered_rnorm_spectrum = function(f, k=2, seed, length_one=FALSE) {
  if (! missing(seed)) set.seed(seed)
  f = apply(f,2,as.integer)
  if (class(f) == "matrix") f = as.data.frame(f)
  size = max(sapply(f,max))
  D = ncol(f)
  ftot = expand_seq(rep(size,D),function(n) -n:n,matrix=FALSE)
  nm = paste0("F",seq_len(D))
  names(f) = nm
  names(ftot) = nm
  fmax = do.call(pmax, lapply(ftot,abs))
  ord = do.call(order,c(list(fmax), ftot))
  ftot = ftot[ord,]
  totsize = nrow(ftot)
  
  plot(ftot,type="n",asp=1)
  text(ftot, labels = seq_len(nrow(ftot)))
  
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
    W[fmerge$j == 1,,drop=FALSE] = 0
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
#' @param spectrum power spectrum of the fields (function of frequency)
#' @param corr.profile correlation profile between (function of wave length)
#' @param closed the probability of fields touching
#' @param gap mean gap between upper and lower fields (overwrites closed)
#' @param seed random seed
#' @param length_one use complex numbers of length one for field generation
#' @param bonds upper and lower bonds for centering
#' @param cut if TRUE, the fields will be cut to bonds
#' @param widen the level of widen effect (see Description)
#' @param widen_grad the inclination of the widen effect
#' 
#' @export
fracture_matrix = function(
    dims = c(10,10),
    span = diag(nrow=length(dims)),
    period = diag(nrow=length(dims)),
    power.spectrum = exp.spectrum(),
    corr.profile = function(k) 0,
    closed = 0.1, gap, seed, length_one = FALSE, bonds, cut=TRUE, widen=0, widen_grad=1) {
  p_ = as.matrix(do.call(expand.grid,lapply(dims,seq_1)))
  p = p_ %*% span
  f_ = as.matrix(do.call(expand.grid,lapply(dims,seq_circ)))
  f = f_ %*% solve(t(span))
  f_len = sqrt(rowSums(f^2))
  
  f_per = fA %*% t(period)
  sel = rowSums(abs(round(f_per) - f_per) > 1e-6) == 0
  coef = matrix(0, nrow(f_per), 2)
  coef[sel,] = ordered_rnorm_spectrum(f_per[sel,,drop=FALSE], k = 2, seed = seed, length_one = length_one)

  Power = power.spectrum(f_len)
  Power[is.infinite(Power)] = 0
    
  coef * Power
  
  Scales = data.frame(
    i  = (1:nrow(fA))[sel],
    kx = as.integer(fA[sel,1]),
    ky = as.integer(fA[sel,2]))
  Scales$k = sqrt(Scales$kx^2 + Scales$ky^2)
  Scales = Scales[Scales$k != 0,]
  Scales$Power = spectrum(Scales$k)
  
  MaxK = max(Scales$kx, Scales$ky, -Scales$kx, -Scales$ky)
  
  RN = ordered_rnorm_mat(MaxK*2+1,MaxK*2+1,3,seed=seed,length_one = length_one)
  
  Scales$RNIndex = ifelse(Scales$kx<0, MaxK*2+1 +Scales$kx, Scales$kx) + ifelse(Scales$ky<0, MaxK*2+1 +Scales$ky, Scales$ky)*(MaxK*2+1) + 1
  Scales$corr = corr.profile(1/Scales$k)
  Scales$corr.angle = atan(Scales$corr)
  Scales$PowerR1F1 = Scales$Power * cos(Scales$corr.angle)
  Scales$PowerR2F1 = Scales$Power * sin(Scales$corr.angle)
  Scales$PowerR1F2 = Scales$Power * sin(Scales$corr.angle)
  Scales$PowerR2F2 = Scales$Power * cos(Scales$corr.angle)
  
  K = matrix(0,N,M)
  K[Scales$i] = Scales$PowerR1F1 * RN[[1]][Scales$RNIndex] + Scales$PowerR2F1 * RN[[2]][Scales$RNIndex]
  f1 = Re(fft(K,inverse=TRUE))
  K[Scales$i] = Scales$PowerR1F2 * RN[[1]][Scales$RNIndex] + Scales$PowerR2F2 * RN[[2]][Scales$RNIndex]
  f2 = Re(fft(K,inverse=TRUE))
  
  Scales$PowerR1Diff = Scales$PowerR1F1 - Scales$PowerR1F2
  Scales$PowerR2Diff = Scales$PowerR2F1 - Scales$PowerR2F2
  diff_sd = sqrt(sum(Scales$PowerR1Diff^2 + Scales$PowerR2Diff^2))
  f1_sd = sqrt(sum(Scales$PowerR1F1^2 + Scales$PowerR2F1^2))
  f2_sd = sqrt(sum(Scales$PowerR1F2^2 + Scales$PowerR2F2^2))
  f1_mean = 0
  f2_mean = 0
  if (is.na(gap)) gap = -qnorm(closed,mean=0,sd=diff_sd)
  f1 = f1 + gap/2
  f1_mean = f1_mean + gap/2
  f2 = f2 - gap/2
  f2_mean = f2_mean - gap/2
  
  sel = f1 < f2
  fm = (f1 + f2)/2
  f1[sel] = fm[sel]
  f2[sel] = fm[sel]
  h2 = f1 - fm
  widen_fac = (pnorm(h2,sd=1/(sqrt(2*pi)*widen_grad/(2*widen)))*2-1)*widen
  f1 = f1 + widen_fac
  f2 = f2 - widen_fac
  if (!missing(bonds)) { if (length(bonds) == 2) {
    shift = mean(bonds)
    f1 = f1 + shift
    f1_mean = f1_mean + shift
    f2 = f2 + shift
    f2_mean = f2_mean + shift
    Pcut = c(pnorm(bonds,mean=f1_mean,sd=f1_sd),pnorm(bonds,mean=f2_mean,sd=f2_sd))
    Pcut = max(Pcut[1],1-Pcut[2],Pcut[3],1-Pcut[4])
    cat("Probability of being out of bonds:", Pcut,"\n")
    if (cut) {
      if (Pcut > 1e-2) warning("Probability of being out of bonds is higher then 1%")
      f1[] = ifelse(f1<bonds[1],bonds[1],ifelse(f1>bonds[2],bonds[2],f1))
      f2[] = ifelse(f2<bonds[1],bonds[1],ifelse(f2>bonds[2],bonds[2],f2))
    }
  } else stop("bonds should be of length 2")}
  ret = list(f1 = f1, f2 = f2, points = pA, scales = Scales, mat=mat, A=A, f1_mean=f1_mean, f1_sd=f1_sd, f2_mean=f2_mean, f2_sd=f2_sd)
  class(ret) = "fracture_field"
  ret
}


