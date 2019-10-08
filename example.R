source("seam_geom.R")
library(rgl)

refine = 64
obj = seam_geom(refine,seed=123)
seam3d(seam.cut(obj))

fname = sprintf("~/seam0_R%02d",refine)

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
