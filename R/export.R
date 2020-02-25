
save.msh = function(obj, filename,type=c("top","bottom"),top="top" %in% type,bottom="bottom" %in% type,middle="middle" %in% type) {
  points = NULL
  triangles = NULL
  i = 0
  n = nrow(obj$points)
  if (top) {
    points = rbind(points, obj$points[,c("f1","x","y")])
    triangles = rbind(triangles, obj$triangles + i*n)
    i = i + 1
  }
  if (bottom) {
    points = rbind(points, obj$points[,c("f2","x","y")])
    triangles = rbind(triangles, obj$triangles + i*n)
    i = i + 1
  }
  if (middle) {
    points = rbind(points, obj$points[,c("fm","x","y")])
    triangles = rbind(triangles, obj$triangles + i*n)
    i = i + 1
  }
  f = file(filename,"w")
  cat("Triangles\n",file=f)
  cat("3D-Nodes ",nrow(obj$points),"\n",sep="",file=f)
  i = 1:nrow(points)-1
  write.table(cbind(i,i,0,points),file=f,row.names=FALSE,col.names=FALSE,sep="\t")
  cat("\nTri3 ",nrow(triangles),"\n",file=f,sep="")
  i = 1:nrow(triangles)-1
  write.table(cbind(i,0,triangles-1),file=f,row.names=FALSE,col.names=FALSE,sep="\t")
  close(f)
}
