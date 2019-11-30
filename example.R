library(rgl)

source("seam_geom.R")

basename = "seam/seam4"
dir.create("seam")
shape=4
s=1
refine=16
for (s in 1:1) {
for (refine_s in c(32)) {
  refine = refine_s * s
  spec = exp.spectrum(scale = 0.015*s^(2.2-2), alpha = 2.2)
  spec_die = ogilve.corr.profile(ML=1.15/s, TL=0.3/s, MaxMF=0, MinMF=1)
  corr = ogilve.corr.profile(ML=0.7/s, TL=0.3/s, MaxMF=1, MinMF=0.1)
  obj = seam.geom(
    refine = refine,
    touch = "include",
    seed=124,
    spectrum = function(k) spec(k)*spec_die(1/k),
    corr.profile = corr,
    closed = 0.1,
    bonds=c(0,0.2)/s,
    widen = 0.01/s,
    widen_grad = 3
  )
  obj = seam.cut(obj)
  obj$points = obj$points*s

  obj_a = seam.touching(obj)
  obj_b = seam.touching(obj,touch = "only")
  fname = sprintf("%s_S%02d_R%02d",basename, s, refine_s)
  seam3d(obj_a)
  writeSTL(paste0(fname,".stl"))
  seam3d(obj_a,type="top")
  writeSTL(paste0(fname,"_top.stl"))
  seam3d(obj_a,type="bottom")
  writeSTL(paste0(fname,"_bottom.stl"))
  seam3d(obj_a,type="middle")
  writeSTL(paste0(fname,"_middle.stl"))
  seam3d(obj_b,type="top")
  writeSTL(paste0(fname,"_touch.stl"))
  border3d(obj)
  writeSTL(paste0(fname,"_border.stl"))
  l = -0.4
  border3d(obj,f1 = l,add=FALSE,col=3)
  triangles3d(t(matrix(c(l,0,0,l,1,0,l,0,1,l,0,1,l,1,0,l,1,1),3,6)),col=3)
  writeSTL(paste0(fname,"_border_b.stl"))
  l = 0.2+0.4
  border3d(obj,f2 = l,add=FALSE,col=3)
  triangles3d(t(matrix(c(l,0,0,l,1,0,l,0,1,l,0,1,l,1,0,l,1,1),3,6)),col=3)
  writeSTL(paste0(fname,"_border_t.stl"))
}
}
source("seam_geom.R")

refine = 32
  #  obj = seam.geom(refine,seed=234,generator=seam.geom.my())
  obj_t = seam.geom(
    refine = refine,
    touch = "only",
    seed=124,
    spectrum=exp.spectrum(scale = 0.015, alpha = 2.2),
    corr.profile = ogilve.corr.profile(ML=0.7, TL=0.3, MaxMF=0.9, MinMF=0.1),
    closed = 0.1,
    bonds=c(0,0.2)
  )
  obj = seam.geom(
    refine = refine,
    touch = "exclude",
    seed=124,
    spectrum=exp.spectrum(scale = 0.015, alpha = 2.2),
    corr.profile = ogilve.corr.profile(ML=0.7, TL=0.3, MaxMF=0.9, MinMF=0.1),
    closed = 0.1,
    bonds=c(0,0.2)
  )
  fname = sprintf("%s_R%02d",basename, refine)
  obj_cut = seam.cut(obj)
  seam3d(obj_cut,type = "top")
  writeSTL(paste0(fname,"_cut_t.stl"))
  seam3d(obj_cut,type = "bottom")
  writeSTL(paste0(fname,"_cut_b.stl"))

  seam3d(obj_cut,type = "bottom")
  seam3d(seam.cut(obj_t),type = "bottom",col=2,add=TRUE)
  #writeSTL(paste0(fname,"_cut_b.stl"))
  writeOBJ(paste0(fname,"_cut_b.obj"))

  
  
  source("seam_geom.R")
  refine = 32*4
  obj = seam.geom(
    refine = refine,
    touch = "exclude",
    seed=124,
    spectrum=exp.spectrum(scale = 0.015, alpha = 2.2),
    corr.profile = ogilve.corr.profile(ML=0.7, TL=0.3, MaxMF=0.9, MinMF=0.1),
    closed = 0.1,
    bonds=c(0,0.2),
    widen = 0.005, widen_grad = 2
  )
  seam3d(seam.cut(obj),type = "bottom")  
  obj_t = seam.geom(
    refine = refine,
    touch = "only",
    seed=124,
    spectrum=exp.spectrum(scale = 0.015, alpha = 2.2),
    corr.profile = ogilve.corr.profile(ML=0.7, TL=0.3, MaxMF=0.9, MinMF=0.1),
    closed = 0.1,
    bonds=c(0,0.2),
    widen = 0.005, widen_grad = 3
  )
  seam3d(seam.cut(obj_t),type = "bottom",add = TRUE,col=2)
  fname = sprintf("%s_R%02d",basename, refine)
  writeOBJ(paste0(fname,"_cut_b.obj"))
  
  
  seam3d(seam.cut(obj))  
  
  
  
  refine = 32
  obj = seam.geom(
    refine = refine,
    touch = "exclude",
    seed=124,
    spectrum = exp.spectrum(scale = 0.015, alpha = 2.2),
    corr.profile = ogilve.corr.profile(ML=0.7, TL=0.3, MaxMF=0.9, MinMF=0.1),
    closed = 0.1,
    bonds=c(0,0.2),
    widen = 0.005, widen_grad = 2
  )
  seam3d(seam.cut(obj))  
  

for (s in 1:4) {
  refine = 16*s
  spec = exp.spectrum(scale = 0.015*s^(2.2-2), alpha = 2.2)
  spec_die = ogilve.corr.profile(ML=1.15/s, TL=0.3/s, MaxMF=0, MinMF=1)
  corr = ogilve.corr.profile(ML=0.7/s, TL=0.3/s, MaxMF=1, MinMF=0.1)
  obj = seam.geom(
    refine = refine,
    touch = "exclude",
    seed=121,
    spectrum = function(k) { spec(k) * spec_die(1/k) },
    corr.profile = corr,
    closed = 0.1#, bonds=c(0.0,0.2)/s
  )
  obj = seam.cut(obj)
  obj$points = obj$points*s
  
  seam3d(obj)  
  writeSTL(sprintf("test%02d.stl",s))
  print(range(obj$points$f1,obj$points$f2))
  print(mean(obj$points$h))
}
  
  
refine_s = 32
s = 1
  refine = refine_s * s
  spec = exp.spectrum(scale = 0.015*s^(2.2-2), alpha = 2.2)
  spec_die = ogilve.corr.profile(ML=1.15/s, TL=0.3/s, MaxMF=0, MinMF=1)
  corr = ogilve.corr.profile(ML=0.7/s, TL=0.3/s, MaxMF=1, MinMF=0.1)
  obj = seam.geom(
    refine = refine,
    touch = "include",
    seed=121,
    spectrum = function(k) spec(k)*spec_die(1/k),
    corr.profile = corr,
    closed = 0.1,
    bonds=c(0,0.2)/s,
    widen = 0.01/s,
    widen_grad = 3
  )
  obj = seam.cut(obj)
  obj$points = obj$points*s
  
  fname = sprintf("%s_S%02d_R%02d",basename, s, refine_s)
  seam3d(obj)
  
  border3d(obj,f1 = 1, col=4)

  
  
        
  source("hip.R")
  source("seam_balls.R")
Rmin = 0.004
Rmax = 0.02

dist = function(k,Rmin,Rmax) { rhip(k,Rmin^3,Rmax^3)^(1/3) }
plot(dist(1000,Rmin,Rmax))

margin=0.0005
Rmin = 0.005
Rmax = 0.012
B_L = seam.balls(obj_a,6000,Rmax = Rmax,Rmin=Rmin,dist=dist, mean.neighbor = 500,iterations = 40,margin = margin,max_add = 100)
B = B_L

Rmin = 0.005
Rmax = 0.005
B = seam.balls(obj_a,B=B_L,1000,Rmax = Rmax,Rmin=Rmin,dist=dist, mean.neighbor = 500,iterations = 40,margin = margin,max_add = 500)


sum(4/2*pi*B$r^3)/ seam.volume(obj)

sel = B$x > B$f1 | B$x < B$f2
sum(sel)
B = B[!sel,]

seam3d(seam.cut(obj))
spheres3d(B$x, B$y, B$z, col = "green", radius = B$r)
#triangles3d(P$f2[iv], P$x[iv], P$y[iv], col = 1,alpha=0.7)

clear3d()
spheres3d(B$x, B$y, B$z, col = "green", radius = B$r)


write.csv(B[,c("x","y","z","r")],file=paste0(fname,"_balls.csv"),row.names=FALSE)
write.lammps.data(B,filename=paste0(fname,"_balls.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))

x = seq(0,1,len=nrow(B))
plot(x, sort(B$r,decreasing = TRUE)^3)
points(x, qhip(x, Rmax^3, Rmin^3),col=3)

x = seq(0,1,len=nrow(B))
plot(x, sort(B$r,decreasing = TRUE))
points(x, qhip(x, Rmax, Rmin),col=3)


sum(4/3*pi*B$r^3)/seam.volume(seam.cut(obj))

range(obj$points$f1,obj$points$f2)

source("seam_geom.R")

obj = seam.geom(
  refine = refine,
  seed=234,
  generator = seam.geom.brown(
    scale=0.04,
    alpha = 2.3,
    corr.profile = ogilve.corr.profile(ML=0.5, TL=0.3, MaxMF=0.9, MinMF=0.1),
    closed = 0.1,
    beta.make=TRUE
  )
)
seam3d(seam.cut(obj))




obj = seam.geom(
  refine = refine,
  seed=234,
  generator = seam.geom.brown(
    scale=0.024,
    alpha = 2.3,
    corr.profile = ogilve.corr.profile(ML=0.7, TL=0.3, MaxMF=0.9, MinMF=0.1),
    closed = 0.1,
    beta.make=TRUE,
    beta.shape = 3
  )
)
seam3d(seam.cut(obj))


obj = seam.geom(
  refine = refine,
  seed=234,
  generator = seam.geom.brown(
    scale=0.04,
    alpha = 2.3,
    corr.profile = ogilve.corr.profile(ML=0.5, TL=0.3, MaxMF=0.9, MinMF=0.1),
    closed = 0.1,
    beta.make=FALSE
  )
)
seam3d(seam.cut(obj))

