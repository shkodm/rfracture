#' Expotential spectrum
#' 
#' @param scale The scale factor
#' @param alpha The exponent
#' @param length The reference wavelength
#' 
#' @export
exp.spectrum = function(scale=1, alpha=2, length=1) {
  function(f) scale^2/((f*length)^alpha)
}
