library(rgl)

source("seam_geom.R")

basename = "~/seam/seam4"
shape=4
s=1
refine_s=32
for (s in 1:1) {
for (refine_s in c(16,64,32)) {
  refine = refine_s * s
  spec = exp.spectrum(scale = 0.015*s^(2.2-2), alpha = 2.2)
  spec_die = ogilve.corr.profile(ML=0.95/s, TL=0.1/s, MaxMF=0, MinMF=1)
  corr = ogilve.corr.profile(ML=0.7/s, TL=0.3/s, MaxMF=1, MinMF=0.1)
  obj = seam.geom(
    refine = refine,
    seed=124,
    spectrum = function(k) spec(k)*spec_die(1/k),
    corr.profile = corr,
    closed = 0.1,
    bonds=c(0,0.2)/s, cut= TRUE,
    widen = 0.01/s,
    widen_grad = 3
  )
  obj = seam.cut(obj)
  obj$points = obj$points*s
  print(range(obj$points$f1,obj$points$f2))
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


basename = "~/seam/seamf"
refine_s = 1
s=1
k = s*refine_s
fname = sprintf("%s_S%02d_R%02d",basename, s, refine_s)
I = matrix(1:((k+1)*(k+1)),k+1,k+1)
points = expand.grid(x = (0:k)/refine_s,y = (0:k)/refine_s)
points$f2 = 0
points$f1 = 0.2*31/32
points$h = points$f1 - points$f2
triangles = rbind(
  cbind(as.vector(I[1:k,1:k]),as.vector(I[1:k,1:k+1]),as.vector(I[1:k+1,1:k])),
  cbind(as.vector(I[1:k+1,1:k+1]),as.vector(I[1:k,1:k+1]),as.vector(I[1:k+1,1:k]))
)
obj_a = list(points=points, triangles=triangles)
seam3d(obj_a)
writeSTL(paste0(fname,".stl"))
seam3d(obj_a,type="top")
writeSTL(paste0(fname,"_top.stl"))
seam3d(obj_a,type="bottom")
writeSTL(paste0(fname,"_bottom.stl"))
seam3d(obj_a,type="middle")
writeSTL(paste0(fname,"_middle.stl"))


  source("hip.R")
  source("seam_balls.R")

N = 62500
(seam.volume(seam.cut(obj))*0.2 / (4/3*pi*N))^(1/3)

Rmin = 0.0030
Rmax = 0.0032

seam.volume(seam.cut(obj_a))*0.2 / (4/3*pi*Rmax^3)


dist = function(k,Rmin,Rmax) { rhip(k,Rmin^3,Rmax^3)^(1/3) }
#plot(dist(1000,Rmin,Rmax))

margin=0.0017
B = seam.balls(obj_a,N,Rmax = Rmax,Rmin=Rmin,dist=dist, mean.neighbor = 500,iterations = 40,margin = margin,max_add = 8000,period=s)
#B2 = seam.balls(obj,3000,Rmax = Rmax,Rmin=Rmin,B=B,mean.neighbor = 500,iterations = 40,margin = margin,max_add = 30,seed=7)
#B = B2
BS = B
sel = B$x > B$f1 | B$x < B$f2
sum(sel)
B = B[!sel,]

seam3d(obj_a)
spheres3d(B$x, B$y, B$z, col = "green", radius = B$r)
#triangles3d(P$f2[iv], P$x[iv], P$y[iv], col = 1,alpha=0.7)

seam.volume(seam.cut(obj))*0.2 / (4/3*pi*Rmax^3)
seam.volume(obj)*0.2 / (4/3*pi*Rmax^3)


write.csv(B[,c("x","y","z","r")],file=paste0(fname,"_balls_L.csv"),row.names=FALSE)
write.lammps.data(B,filename=paste0(fname,"_balls_L.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))


B = seam.balls(obj_a,N, B=BS, Rmax = Rmax,Rmin=Rmin,dist=dist, mean.neighbor = 500,iterations = 40,margin = margin,max_add = 8000,period=s)

write.csv(B[,c("x","y","z","r")],file=paste0(fname,"_balls_L2.csv"),row.names=FALSE)
write.lammps.data(B,filename=paste0(fname,"_balls_L2.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))



################

margin=0.0017
N=64000
B = seam.balls(obj_a,N,Rmax = Rmax,Rmin=Rmin,dist=dist, mean.neighbor = 500,iterations = 40, relax_iterations = 0, delete=FALSE, margin = margin,max_add = N,period=s)

write.csv(B[,c("x","y","z","r")],file=paste0(fname,"_balls_O.csv"),row.names=FALSE)
write.lammps.data(B,filename=paste0(fname,"_balls_O.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))

#################

N = 1500
margin=0.0017
Rmin = 0.018
Rmax = 0.023
dist = function(k,Rmin,Rmax) { rhip(k,Rmin^3,Rmax^3)^(1/3) }

B = seam.balls(obj_a, N, Rmax=Rmax, Rmin=Rmin, dist=dist, mean.neighbor=500, iterations=40, relax_iterations=0, delete=FALSE, margin=margin, max_add=N, period=s)

write.csv(B[,c("x","y","z","r")],file=paste0(fname,"_balls_O.csv"),row.names=FALSE)
write.lammps.data(B,filename=paste0(fname,"_balls_O.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))


###############
ret = read.lammps.data("~/seam/pack1/out1.dat")
B=ret$B
seam3d(obj_a)
spheres3d(B$x, B$y, B$z, col = "green", radius = B$r)
write.lammps.data(B,filename=paste0(fname,"_balls_R.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))

rbd = function(dx=0,dy=0,dz=0) {
  B = ret$B
  B$x = B$x + dx
  B$y = B$y + dy
  B$z = B$z + dz
  B
}

set = list(
  expand.grid(dx=0  ,dy=0  ,dz=0  ), # 1
  expand.grid(dx=0  ,dy=0:1,dz=0  ), # 2
  expand.grid(dx=0  ,dy=0:1,dz=0:1), # 4
  expand.grid(dx=0  ,dy=0:2,dz=0:1), # 6
  expand.grid(dx=0  ,dy=0:3,dz=0:1), # 8
  expand.grid(dx=0  ,dy=0:3,dz=0:2), # 12
  expand.grid(dx=0  ,dy=0:3,dz=0:3), # 16
  expand.grid(dx=0  ,dy=0:3,dz=0:4), # 20
  expand.grid(dx=0  ,dy=0:3,dz=0:5), # 24
  expand.grid(dx=0  ,dy=0:5,dz=0:4)  # 30
)


for (mul in set[8:10]) {
  sz = sapply(mul,max)+1
  nm = paste0(sz,collapse="")
  B = NULL
  for (i in 1:nrow(mul)) {
    NB = ret$B
    NB$x = NB$x + mul$dx[i]
    NB$y = NB$y + mul$dy[i]
    NB$z = NB$z + mul$dz[i]
    B = rbind(B, NB)
    obj_d = obj_a
    obj_d$points$x = obj_d$points$x + mul$dy[i]
    obj_d$points$y = obj_d$points$y + mul$dz[i]
    seam3d(obj_d,add = i > 1)
  }
  writeSTL(paste0(fname,"_",nm,".stl"))
  write.lammps.data(B,filename=paste0(fname,"_",nm,"_balls.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,sz[2]),zlim=c(0,sz[3]))
}

B = rbind(rbd(0,0,0))
write.lammps.data(B,filename=paste0(fname,"_balls_R_11.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))
B = rbind(rbd(0,0,0),rbd(0,0,0))
write.lammps.data(B,filename=paste0(fname,"_balls_R_11.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))
B = rbind(rbd(0,0,0))
write.lammps.data(B,filename=paste0(fname,"_balls_R_11.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))
B = rbind(rbd(0,0,0))
write.lammps.data(B,filename=paste0(fname,"_balls_R_11.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))
B = rbind(rbd(0,0,0))
write.lammps.data(B,filename=paste0(fname,"_balls_R_11.dat"),density = 2,xlim = c(0,0.2),ylim=c(0,1),zlim=c(0,1))

4/3*pi*sum(B$r^3)/seam.volume(seam.cut(obj))
diff(obj$points$y[obj$points$x == 1])


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

