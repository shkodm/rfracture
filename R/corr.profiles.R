#' Correlation profile with linear cut-off
#' 
#' @param ML Length of correlation cut-off
#' @param TL Width of linear cut-off fade
#' @param MinMF Minimal correlation
#' @param MaxMF Maximal correlation
#' @references
#' Steven R. Ogilvie, Evgeny Isakov, Paul W.J. Glover (2006), Fluid flow through rough fractures in rocks. II: A new matching model for rough rock fractures,
#' Earth and Planetary Science Letters, Volume 241, Issues 3–4, https://doi.org/10.1016/j.epsl.2005.11.041.
#' 
#' @export
ogilvie.corr.profile = function(ML=0.5, TL=0, MinMF=0, MaxMF=1) function(lambda) {
  ifelse(lambda >= ML+TL/2, 1, ifelse(lambda <= ML-TL/2, 0, (ML-TL/2-lambda)*(lambda-(ML+3*TL/2))/(TL*TL) ))*(MaxMF-MinMF)+MinMF
}
