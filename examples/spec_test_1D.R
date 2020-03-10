library(rfracture)
library(rgl)
library(progress)
library(parallel)

cores=detectCores()
cl <- makeCluster(cores-1)

nx = seq(0,1,len=300)
refine = 2^(0:5)
freq = NULL
method = "diagonals"
power.spectrum = exp_spectrum(scale=0.02,alpha=2)
power.spectrum = function(f) ifelse(f<5, 0, 0.004/(f^2))
#power.spectrum = function(f) 0.02^2*exp(-2*(f/5)^2)
repetitions = 200

tab = expand.grid(refine=refine, gauss.order = 1:3, stringsAsFactors = FALSE)
pb <- progress_bar$new(total = nrow(tab))

sp = lapply(seq_len(nrow(tab)), function(i) {
  refine = tab$refine[i]
  gauss.order = tab$gauss.order[i]
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum, seed=123, method=method)
  pb$tick()
  clusterExport(cl = cl, varlist = c("refine","method","power.spectrum","nx","gauss.order"),envir = environment())
  ny = parSapplyLB(cl, seq_len(repetitions), function(j) {
    library(rfracture)
    ret = fracture_matrix(5*refine, corr.profile=function(lambda) 1,gap=0.05, power.iso=power.spectrum,gauss.order = gauss.order)
    x = c(as.vector(ret$points),1)
    y = c(as.vector(ret$f1),ret$f1[1])
    y
  })
  
  tny = ts(ny, deltat = 1/(5*refine))
  
  sp = spectrum(tny,plot=FALSE)
  sp$spec = rowMeans(sp$spec)
  sp
})
stopCluster(cl)

tab$refine_f = factor(tab$refine)

freq = range(sapply(sp, function(x) range(x$freq)))
speclim = range({x = sapply(sp, function(x) range(x$spec)); x[x<1e-16] = NA; x}, na.rm = TRUE)
freq = seq(freq[1],freq[2],len=200)
plot(freq, power.spectrum(freq),type="l",log="xy",ylim=speclim)
for (i in 1:length(sp)) {
  lines(sp[[i]]$freq, sp[[i]]$spec, col=as.integer(tab$refine_f[i]), pch=tab$gauss.order[i])
}

