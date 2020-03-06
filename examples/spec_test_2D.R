library(rfracture)
library(rgl)
library(progress)
library(parallel)

cores=detectCores()

nx = seq(0,1,len=300)
refine = 2^(0:5)
freq = NULL
method = "diagonals"
power.spectrum = exp_spectrum(scale=0.02,alpha=3.5)
#power.spectrum = function(f) 0.02^2*exp(-2*(f/5)^2)
repetitions = 200

tab = expand.grid(refine=refine, length_one = c(TRUE,FALSE), stringsAsFactors = FALSE)
pb <- progress_bar$new(total = nrow(tab))

sp = lapply(seq_len(nrow(tab)), function(i) {
  refine = tab$refine[i]
  length_one = tab$length_one[i]
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum, seed=123, method=method)
  pb$tick()
  cl <- makeCluster(cores-1)
  clusterExport(cl, c("refine","method","power.spectrum","nx","length_one"), envir = environment())
  spl = parLapplyLB(cl, seq_len(repetitions), function(j) {
    library(rfracture)
    #ret = fracture_matrix(5*refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum,length_one = length_one)
    ret = fracture_matrix(c(1,1)*5*refine, corr.profile=function(lambda) 1, gap=0.05, power.spectrum=power.spectrum, length_one = length_one)
    #x = c(as.vector(ret$points[1:(5*refine),1]),1)
    #y = c(as.vector(ret$f1[,1]),ret$f1[1,1])
    y = ret$f1
    tny = ts(y, deltat = 1/(5*refine))
    sp = spectrum(tny,plot=FALSE)
    sp$spec = rowMeans(sp$spec)
    sp
  })
  stopCluster(cl)
  sp = spl[[1]]
  sp$spec = rowMeans(sapply(spl, function(x) x$spec))
  sp
})

tab$refine_f = factor(tab$refine)
tab$length_one_f = factor(tab$length_one)

freq = range(sapply(sp, function(x) range(x$freq)))
freq = seq(freq[1],freq[2],len=200)
plot(freq, 2*freq*power.spectrum(freq), type="l",log="xy")
for (i in 1:length(sp)) {
  freq = sp[[i]]$freq
  spec = sp[[i]]$spec
  points(freq, spec, col=as.integer(tab$refine_f[i]), pch=as.integer(tab$length_one_f[i]))
  lines(freq, 2*power.spectrum(freq)*freq * (1-(freq/max(freq*1.5))^2),lty=2,col=as.integer(tab$refine_f[i]))
  mfreq = (5*tab$refine[i]/2)
  freq = -mfreq:mfreq
  M = outer(freq, freq, function(x,y) power.spectrum(sqrt(x*x+y*y)))
  lines(freq, colSums(M),col=as.integer(tab$refine_f[i]))
}

legend("bottomleft",legend = c(levels(tab$refine_f), levels(tab$length_one_f))
       ,pch=c(rep(1,nlevels(tab$refine_f)), seq_len(nlevels(tab$length_one_f)))
       ,col=c(seq_len(nlevels(tab$refine_f)), rep(1,nlevels(tab$length_one_f)))
)

for (i in 1:nrow(tab)) {
  mfreq = (5*tab$refine[i]/2)
  freq = -mfreq:mfreq
  M = outer(freq, freq, function(x,y) power.spectrum(sqrt(x*x+y*y)))
  lines(freq, power.spectrum(freq)*2*freq, col=4,lwd=3)
  lines(freq, colSums(M),col=3,lwd=3)
}



refine = 10
ret1 = fracture_matrix(c(1)*5*refine, corr.profile=function(lambda) 1, gap=0.05, power.spectrum=power.spectrum)
ret2 = fracture_matrix(c(1,1)*5*refine, corr.profile=function(lambda) 1, gap=0.05, power.spectrum=power.spectrum)
ret3 = fracture_matrix(c(1,1,1)*5*refine, corr.profile=function(lambda) 1, gap=0.05, power.spectrum=power.spectrum)

matplot(cbind(ret1$f1,ret2$f1[1,]/2,ret3$f1[,1,1]),type="l",lty=1)

c(ret1$var.midline,ret2$var.midline,ret3$var.midline) / ret1$var.midline
sqrt(2*pi)


c(sd(ret1$f1),sd(ret2$f1),sd(ret3$f1)) / sd(ret1$f1)
