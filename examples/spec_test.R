library(rfracture)
library(rgl)
library(progress)
library(parallel)

myspectrum = function(x,... ) {
  mx = as.matrix(x)
  n = nrow(mx)
  span = n/frequency(x)
  list(
    spec = apply(mx,2,function(x) Mod(fft(x)/n)^2)*span,
    freq = abs(seq_circ(n)/span)
  )
}

cores=detectCores()
cl <- makeCluster(cores-1)

nx = seq(-0.25,0.25,len=400)
refine = 2^(0:4)
freq = NULL
method = "diagonals"
#power.spectrum = exp_spectrum(scale=0.02,alpha=3.5)
#power.spectrum = function(f) ifelse(f<10, 0, 0.004/(f^4.5))
#power.spectrum = function(f) 0.004/(f^4.5)
power.spectrum = function(f) 0.004/(f^4)
#power.spectrum = function(f) 0.02^2*exp(-2*(f/5)^2)
repetitions = 80

tab = expand.grid(method=c("triangles", "diagonals"), refine=refine, stringsAsFactors = FALSE)
pb <- progress_bar$new(total = nrow(tab))

sp = lapply(seq_len(nrow(tab)), function(i) {
  refine = tab$refine[i]
  method = tab$method[i]
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=exp_spectrum(scale=0.02,alpha=2.5), seed=123)
  #ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.spectrum=power.spectrum, seed=123, method=method)
  nx = seq(-0.25,0.25,len=refine*5)
  
  clusterExport(cl, c("refine","method","power.spectrum","nx","myspectrum"), envir = environment())
  spl = parLapplyLB(cl, seq_len(repetitions), function(j) {
    library(rfracture)
    
    #ny = replicate(repetitions, {
    ret = fracture_geom(width=1, refine=refine, corr.profile=function(lambda) 1,gap=0.05, power.iso=power.spectrum, method=method, gauss.order = 1, seed=12*j)
    
    d = rnorm(2)
    #d = c(0,rnorm(1))
    d = d/sqrt(sum(d^2))
    ret$points$zeta2 = (ret$points$x - 0.5)*d[1] + (ret$points$y - 0.5)*d[2]
    ret$points$zeta1 = (ret$points$x - 0.5)*d[2] - (ret$points$y - 0.5)*d[1]
    
    ret2 = slice(ret, "zeta2", level = runif(1, -0.25,0.25), value="edge")
    
    x = ret2$points$zeta1
    y = ret2$points$f1
    
    ny = approx(x, y, xout = nx)
    tny = ts(ny$y, deltat = nx[2]-nx[1])
    sp = spectrum(tny,plot=FALSE)
    #sp = myspectrum(tny,plot=FALSE)
    sp
  })
  
  sp = spl[[1]]
  sp$spec = rowMeans(sapply(spl, function(x) x$spec))
  pb$tick()
  sp
})
stopCluster(cl)


tab$refine_f = factor(tab$refine*5)
tab$method_f = factor(tab$method)

freq = range(sapply(sp, function(x) range(x$freq)))
freq = seq(freq[1],freq[2],len=200)
aspec = pi*freq*power.spectrum(freq)
speclim = range(sapply(sp, function(x) range(x$spec)))

pdf("spec2Dlines.pdf")
plot(freq, aspec, type="n", log="y", ylim=speclim, lty=2, xlim=c(0,80), xlab="Frequency [1/m]", ylab="PSD [m2]", yaxt='n')
for (i in 1:length(sp)) {
  freq = sp[[i]]$freq
  spec = sp[[i]]$spec
  points(freq, spec, col=as.integer(tab$refine_f[i]), pch=as.integer(tab$method_f[i]))
}
freq = 1:80
M = outer(freq, freq, function(x,y) power.spectrum(sqrt(x*x+y*y)))*2
lines(freq, colSums(M), lty=2, lwd=2)
a = log10(axTicks(2))
a = floor(min(a)):ceiling(max(a))
axis(2, at=10^a, labels = sapply(a, function(a) as.expression(bquote(10**.(a)))))

legend("topright",legend = c(levels(tab$refine_f), levels(tab$method_f))
       ,pch=c(rep(1,nlevels(tab$refine_f)), seq_len(nlevels(tab$refine_f)))
       ,col=c(seq_len(nlevels(tab$refine_f)), rep(1,nlevels(tab$refine_f)))
       ,title="Resolution"
)
dev.off()
