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

sp = sapply(seq_len(nrow(tab)), function(i) {
  refine = tab$refine[i]
  length_one = tab$length_one[i]
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum, seed=123, method=method)
  pb$tick()
  clusterExport(cl, c("refine","method","power.spectrum","nx"))
  ny = parSapplyLB(cl, seq_len(repetitions), function(j) {
    library(rfracture)
    
    #ny = replicate(repetitions, {
#    ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum, method=method)
    ret = fracture_matrix(5*refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum,length_one = length_one)
    x = c(as.vector(ret$points),1)
    y = c(as.vector(ret$f1),ret$f1[1])
    
    ny = approx(x, y, xout = nx)
    ny$y
  })
  
  tny = ts(ny, nx[1],deltat = nx[2]-nx[1])
  
  sp = spectrum(tny,plot=FALSE)
  
  freq <<- sp$freq
  rowMeans(sp$spec)
})

stopCluster(cl)

tab$refine_f = factor(tab$refine)


matplot(freq, sp, log="xy", col=as.integer(tab$refine_f),type="b", ylim=c(1e-10,0.02^2))

matplot(freq, sp, log="xy", col=as.integer(tab$refine_f),type="b", ylim=c(1e-10,0.02^2))
lines(freq, power.spectrum(freq))
legend("bottomleft",legend = c(levels(tab$refine_f), levels(tab$method_f))
       ,pch=c(rep(1,nlevels(tab$refine_f)), seq_len(nlevels(tab$refine_f)))
       ,col=c(seq_len(nlevels(tab$refine_f)), rep(1,nlevels(tab$refine_f)))
)

