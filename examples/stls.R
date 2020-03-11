
G = 12
mar = 1/G/32 * 0.6
refine = 10
closed = 0.1
for (closed in seq(0,0.5,0.01)[-1]){
  ret = fracture_geom(refine=refine, closed = closed, corr.profile = corr.profile, power.iso = power.iso, seed=125)
  ret2 = slice(ret,  value="above")
  dh = -mean(range(ret2$points$f1, ret2$points$f2)) + 1/(G*2)
  ret$points$f1 = ret$points$f1 + dh
  ret$points$f2 = ret$points$f2 + dh
  ret$points$fm = ret$points$fm + dh
  ret2 = ret
  ret2 = slice(ret2, by = "f1", level=1/G - mar, flatten="above")
  ret2 = slice(ret2, by = "f2", level=  0 + mar, flatten="below")
  ret2 = slice(ret2,  value="above")
  fracture3d(cut(ret2))
  cbind(range(ret2$points$f1, ret2$points$f2), c(mar,1/G-mar), c(0,1/G))
  
  writeSTL(sprintf("seam_test_G%02d_R%03d.stl",round(closed*100),refine))
}
