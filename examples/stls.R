library(rfracture)

alpha = 4.5
power.iso = function(f) 0.000001*ifelse(f<5,1,(f/5)^-alpha)
corr.profile = function(lambda) ifelse(lambda<0.5,0.1,0.9)

G = 12
mar = 1/G/32 * 0.6


closed = 0.05
refine = 20
for (refine in 2^{1:6}){
  ret = fracture_geom(refine=refine, corr.profile = corr.profile, power.iso = power.iso, seed=123)
  for (closed in 0.1){
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

closed = 0.01
refine = 32
for (refine in 2^{1:6}){
  ret = fracture_geom(refine=refine, corr.profile = corr.profile, power.iso = power.iso, seed=123)
  for (closed in seq(0.05,0.95,0.05)){
    ret2 = set_gap(ret, closed=closed)
    ret2 = slice(ret2,  value="above")
    ret2 = slice(ret2,  level=1/G - 2*mar, flatten="above")
    if (require(rgl)) fracture3d(cut(ret2),edge.col = 4, vertex.col = 4)
    print(cbind(range(ret2$points$f1, ret2$points$f2), c(mar,1/G-mar), c(0,1/G)))
    write.stl(cut(ret2), sprintf("seam_test_G%02d_R%02d.stl",round(closed*100),refine))
  }
}
ret2$gap*G*32






alpha = 4.5
sd = 0.001
#sd = 0.000
power.iso = function(f) sd^2*ifelse(f<5,1,(f/5)^-alpha)
corr.profile = function(lambda) ifelse(lambda<0.5,0.1,0.9)
corr.profile = function(lambda) 0.2

G = 12
mar = 1/G/32 * 0.6

refine = 32
gap_n = 16
for (refine in 2^{1:6}){
  ret = fracture_geom(refine=refine, corr.profile = corr.profile, power.iso = power.iso, seed=123)
  for (gap_n in 0:16){
    gap = gap_n / (G*32)
    ret2 = set_gap(ret, gap = gap)
    ret2 = slice(ret2,  value="above")
    ret2 = slice(ret2,  level=1/G - 2*mar, flatten="above")
    
    if (require(rgl)) fracture3d(cut(ret2),edge.col = 4, vertex.col = 4)
    print(cbind(range(ret2$points$f1, ret2$points$f2), c(mar,1/G-mar), c(0,1/G)))
    write.stl(cut(ret2), sprintf("seam_test_G%02d_R%02d.stl",gap_n,refine))
  }
}

cov2cor(ret$cov.final)




alpha = 4.5
sd = 0.001
#sd = 0.000
power.iso = function(f) sd^2*ifelse(f<5,1,(f/5)^-alpha)
corr.profile = function(lambda) ifelse(lambda<0.5,0.1,0.9)
corr.profile = function(lambda) 0.9

G = 12
mar = 1/G/32 * 0.6

refine = 32
gap_n = 16
ret = fracture_geom(refine=refine, corr.profile = corr.profile, power.iso = power.iso, seed=123)
ret$points$f1 =  -cos((ret$points$x)*2*pi)*0.025
ret$points$f2 = -ret$points$f1
ret$points$fm = (ret$points$f1 + ret$points$f2)/2
ret$points$h  =  ret$points$f1 - ret$points$f2
gap = gap_n / (G*32)
ret2 = set_gap(ret, gap = gap)
ret2 = slice(ret2,  value="above")
ret2 = cut(ret2)
if (require(rgl)) fracture3d(ret2,edge.col = 4, vertex.col = 4)

res = solve_reynolds(ret2)
res

eigen(res$perm)
