library(rfracture)
library(rgl)
library(progress)
library(parallel)

cores=detectCores()
cl <- makeCluster(cores-1)

nx = seq(-0.25,0.25,len=300)
refine = 2^(0:4)
freq = NULL
method = "diagonals"
power.spectrum = exp_spectrum(scale=0.02,alpha=3.5)
#power.spectrum = function(f) 0.02^2*exp(-2*(f/5)^2)
repetitions = 200

tab = expand.grid(method=c("triangles", "diagonals"), refine=refine, stringsAsFactors = FALSE)
pb <- progress_bar$new(total = nrow(tab))

sp = sapply(seq_len(nrow(tab)), function(i) {
  refine = tab$refine[i]
  method = tab$method[i]
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=exp_spectrum(scale=0.02,alpha=2.5), seed=123)
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum, seed=123, method=method)
  pb$tick()
  clusterExport(cl, c("refine","method","power.spectrum","nx"), envir = environment())
  ny = parSapplyLB(cl, seq_len(repetitions), function(j) {
    library(rfracture)
    
    #ny = replicate(repetitions, {
    ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.iso=power.spectrum, method=method, gauss.order = 3)
    
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
  
  tny = ts(ny, nx[1],deltat = nx[2]-nx[1])
  
  sp = spectrum(tny,plot=FALSE)
  
  freq <<- sp$freq
  rowMeans(sp$spec)
})

stopCluster(cl)

tab$refine_f = factor(tab$refine)
tab$method_f = factor(tab$method)
matplot(freq, sp, log="xy", col=as.integer(tab$refine_f), pch=as.integer(tab$method_f), type="b")

lines(freq, power.spectrum(freq)*freq/pi)
legend("bottomleft",legend = c(levels(tab$refine_f), levels(tab$method_f))
       ,pch=c(rep(1,nlevels(tab$refine_f)), seq_len(nlevels(tab$refine_f)))
       ,col=c(seq_len(nlevels(tab$refine_f)), rep(1,nlevels(tab$refine_f)))
)

abline(v=5*refine/2)

n = 10
freq = 1:n-1
M = outer(freq, freq, function(x,y) power.spectrum(sqrt(x*x+y*y)))
lines(freq, power.spectrum(freq)*freq/pi, col=4,lwd=3)
lines(freq, colSums(M),col=3,lwd=3)

for (i in 1:nrow(tab)) {
  freq = 1:(5*tab$refine[i]/2)-1
  M = outer(freq, freq, function(x,y) power.spectrum(sqrt(x*x+y*y)))
  lines(freq, power.spectrum(freq)*freq/pi, col=4,lwd=3)
  lines(freq, colSums(M),col=3,lwd=3)
}

matplot(freq, sp, log="xy", col=as.integer(tab$refine_f),type="n", ylim=c(1e-10,0.02^2))
for (i in 1:ncol(sp)) {
  #points(freq, sp[,i], pch=ifelse(freq <= 5*tab$refine[i]/2,15,0)+1, col=as.integer(tab$refine_f[i]))
  qs = 5*tab$refine[i]/2
  sel = freq <= qs
  lines(freq[sel], sp[sel,i], col=as.integer(tab$refine_f[i]))
  lines(freq, 2*freq*power.spectrum(freq) * sqrt(1 - (freq/qs)^2))
  #points(freq[sel], sp[sel,i], col=as.integer(tab$refine_f[i]))
}
lines(freq, power.spectrum(freq))





