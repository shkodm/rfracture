library(rgl)

source("seam_geom.R")

basename = "~/seam/seam1"
shape=4
for (refine in c(1,2,4,8,16,32)) {
#  obj = seam.geom(refine,seed=234,generator=seam.geom.my())
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
  seam.volume(seam.cut(obj))
  
  fname = sprintf("%s_R%02d",basename, refine)
  seam3d(obj)
  writeSTL(paste0(fname,".stl"))
  seam3d(obj,type="top")
  writeSTL(paste0(fname,"_top.stl"))
  seam3d(obj,type="bottom")
  writeSTL(paste0(fname,"_bottom.stl"))
  seam3d(obj,type="middle")
  writeSTL(paste0(fname,"_middle.stl"))
  save.msh(obj, paste0(fname,"_top.msh"), type="top")
  save.msh(obj, paste0(fname,"_bottom.msh"), type="bottom")
  seam3d(seam.cut(obj))
  writeSTL(paste0(fname,"_cut.stl"))
}


source("seam_balls.R")
Rmax = 0.02
Rmin = 0.02
B = seam.balls(obj,3000,Rmax = Rmax,Rmin=Rmin,mean.neighbor = 500,iterations = 40,margin = 0,margin_opt = 0,max_add = 30)
#tds = B[[2]]
#X = B[[3]]
#B = B[[1]]
#points3d(X)


seam3d(seam.cut(obj))
spheres3d(B$x, B$y, B$z, col = "green", radius = B$r)
#triangles3d(P$f2[iv], P$x[iv], P$y[iv], col = 1,alpha=0.7)

write.csv(B[,c("x","y","z","r")],file=paste0(fname,"_balls.csv"),row.names=FALSE)
write.lammps.data(B,filename=paste0(fname,"_balls.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))






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

