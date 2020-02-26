#' Generic function for volume calculation
#' 
#' @param x object of which to calculate volume
#' @param ... other parameters
#' 
#' @export
volume  = function (x, ...) UseMethod("volume")
