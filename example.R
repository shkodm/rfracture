source("seam_geom.R")
library(rgl)

shape=4
for (refine in c(1,2,4,8,16,32)) {
  obj = seam.geom(refine,seed=234,shape1=shape,shape2=shape)
  seam3d(seam.cut(obj))
  seam.volume(seam.cut(obj))
  
  fname = sprintf("~/seam/seam0_R%02d",refine)
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

B = seam.balls(obj,3000)
seam3d(seam.cut(obj))
spheres3d(B$x, B$y, B$z, col = "green", radius = B$r)
#triangles3d(P$f2[iv], P$x[iv], P$y[iv], col = 1,alpha=0.7)

write.csv(B[,c("x","y","z","r")],file=paste0(fname,"_balls.csv"),row.names=FALSE)
write.lammps.data(B,filename=paste0(fname,"_balls.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))

