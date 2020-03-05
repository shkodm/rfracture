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
#power.spectrum = function(f) 0.02^2*exp(-2*(f/5)^2)
repetitions = 200

tab = expand.grid(refine=refine, length_one = c(TRUE,FALSE), stringsAsFactors = FALSE)
pb <- progress_bar$new(total = nrow(tab))

sp = lapply(seq_len(nrow(tab)), function(i) {
  refine = tab$refine[i]
  length_one = tab$length_one[i]
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum, seed=123, method=method)
  pb$tick()
  clusterExport(cl, c("refine","method","power.spectrum","nx","length_one"))
  ny = parSapplyLB(cl, seq_len(repetitions), function(j) {
    library(rfracture)
    ret = fracture_matrix(5*refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum,length_one = length_one)
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
freq = seq(freq[1],freq[2],len=200)
plot(freq, power.spectrum(freq),type="l",log="y")
for (i in 1:length(sp)) {
  lines(sp[[i]]$freq, sp[[i]]$spec,col=as.integer(tab$refine_f[i]))
}

