library(rfracture)
library(rgl)
library(progress)
library(parallel)

myspectrum = function(x,... ) {
  n = nrow(x)
  span = n/frequency(x)
  list(
    spec = apply(x,2,function(x) Mod(fft(x)/n)^2)*span,
    freq = seq_circ(n)/span
  )
}

cores=detectCores()
cl <- makeCluster(cores-1)

nx = seq(-0.25,0.25,len=300)
refine = 2^(0:4)
freq = NULL
method = "diagonals"
power.spectrum = exp_spectrum(scale=0.02,alpha=3.5)
power.spectrum = function(f) ifelse(f<10, 0, 0.004/(f^4.5))
power.spectrum = function(f) 0.004/(f^4.5)
#power.spectrum = function(f) 0.02^2*exp(-2*(f/5)^2)
repetitions = 2

tab = expand.grid(method=c("triangles", "diagonals"), refine=refine, stringsAsFactors = FALSE)
pb <- progress_bar$new(total = nrow(tab))

sp = sapply(seq_len(nrow(tab)), function(i) {
  refine = tab$refine[i]
  method = tab$method[i]
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=exp_spectrum(scale=0.02,alpha=2.5), seed=123)
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum, seed=123, method=method)
  clusterExport(cl, c("refine","method","power.spectrum","nx"), envir = environment())
  ny = parSapplyLB(cl, seq_len(repetitions), function(j) {
    library(rfracture)
    
    #ny = replicate(repetitions, {
    ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.iso=power.spectrum, method=method, gauss.order = 1, seed=12*j)
    
    #d = rnorm(2)
    d = c(0,rnorm(1))
    d = d/sqrt(sum(d^2))
    ret$points$zeta2 = (ret$points$x - 0.5)*d[1] + (ret$points$y - 0.5)*d[2]
    ret$points$zeta1 = (ret$points$x - 0.5)*d[2] - (ret$points$y - 0.5)*d[1]
    
    ret2 = slice(ret, "zeta2", level = runif(1, -0.25,0.25),  value="edge")
    
    x = ret2$points$zeta1
    y = ret2$points$f1
    
    ny = approx(x, y, xout = nx)
    ny$y
  })
  
  tny = ts(ny, start = nx[1], deltat = nx[2]-nx[1])
#  sp = spectrum(tny,plot=FALSE)
  sp = myspectrum(tny)
  
  freq <<- sp$freq
  pb$tick()
  rowMeans(sp$spec)
})

stopCluster(cl)

tab$refine_f = factor(tab$refine)
tab$method_f = factor(tab$method)
matplot(abs(freq), sp, log="xy", col=as.integer(tab$refine_f), pch=as.integer(tab$method_f), type="b",lty=1)

lines(freq, power.spectrum(freq)*freq/pi*30)
#lines(freq, power.spectrum(freq))
legend("bottomleft",legend = c(levels(tab$refine_f), levels(tab$method_f))
       ,pch=c(rep(1,nlevels(tab$refine_f)), seq_len(nlevels(tab$refine_f)))
       ,col=c(seq_len(nlevels(tab$refine_f)), rep(1,nlevels(tab$refine_f)))
)


lines(freq, power.spectrum(freq)*freq*freq*3)


abline(v=5*refine/2)

for (i in 1:nrow(tab)) {
  freq_ = 5*tab$refine[i]
  freq_ = (-freq_+1):(freq_)
  M = outer(freq_, freq_, function(x,y) power.spectrum(sqrt(x*x+y*y)))
#  lines(freq, power.spectrum(freq)*freq/pi, col=4,lwd=3)
  lines(freq_, colSums(M),lwd=3, col=as.integer(tab$refine_f)[i])
}


plot(freq, power.spectrum(freq)*freq/pi, log="xy",ylim=c(1e-12,1), type="l")
for (i in 1:nrow(tab)) {
  sel = freq <= 5*tab$refine[i]/2
  points(freq[sel], sp[sel,i], log="xy", col=as.integer(tab$refine_f[i]), pch=as.integer(tab$method_f[i]), type="b")
}

