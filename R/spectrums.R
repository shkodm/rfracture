#' Expotential spectrum
#' 
#' @param scale The scale factor
#' @param fractal.dimension Fractal dimension of the spectrum
#' @param alpha The exponent
#' @param length The reference wavelength
#' 
#' @export
exp_spectrum = function(scale=1, fractal.dimension=2.5, alpha=7-2*fractal.dimension, length=1) {
  function(f) scale^2/((f*length)^alpha)
}
