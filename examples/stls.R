library(rfracture)

alpha = 4.5
power.iso = function(f) 0.000001*ifelse(f<5,0,(f/5)^-alpha)
corr.profile = function(lambda) ifelse(lambda<0.5,0.1,0.9)

G = 12
mar = 1/G/32 * 0.6

closed = 0.1
for (refine in 2^{1:6}){
  ret = fracture_geom(refine=refine, corr.profile = corr.profile, power.iso = power.iso, seed=125)
  for (closed in 0.2){
    ret2 = set_gap(ret, closed=closed)
    ret2 = slice(ret2,  value="above")
    dh = -mean(range(ret2$points$f1, ret2$points$f2)) + 1/(G*2)
    ret2 = set_gap(ret2, offset=dh)
    ret2 = slice(ret2, by = "f1", level=1/G - mar, flatten="above")
    ret2 = slice(ret2, by = "f2", level=  0 + mar, flatten="below")
    if (require(rgl)) fracture3d(cut(ret2),edge.col = 4, vertex.col = 4)
    print(cbind(range(ret2$points$f1, ret2$points$f2), c(mar,1/G-mar), c(0,1/G)))
    write.stl(cut(ret2), sprintf("seam_test_G%02d_R%02d.stl",round(closed*100),refine))
  }
}

