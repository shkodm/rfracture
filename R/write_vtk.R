obj = ret

write.vtk.tet.mesh = function(obj, filename){
  points = rbind(
    as.matrix(obj$points[,c("f1","x","y")]),
    as.matrix(obj$points[,c("f2","x","y")])
  )
  
  k = nrow(obj$points)
  tr = obj$triangles
  tr = matrix(tr[order(row(tr), tr)], ncol=ncol(tr), byrow = TRUE)
  
  tets = cbind( c(tr[,1]+k,tr[,2]+k,tr[,3]+k),
                c(tr[,1]  ,tr[,2]  ,tr[,3]  ),
                c(tr[,2]  ,tr[,3]  ,tr[,1]+k),
                c(tr[,3]  ,tr[,1]+k,tr[,2]+k))
  
  sel = abs(points[tets[,1],1] - points[tets[,2],1]) < 1e-10
  tets = tets[!sel,]
  
  f = file(filename,"w")
  writeLines(c(
    "# vtk DataFile Version 2.0",
    "Unstructured Grid",
    "ASCII",
    "DATASET UNSTRUCTURED_GRID",
    paste("POINTS",nrow(points),"double")),con=f)
  write.table(points, file=f, row.names = FALSE, col.names=FALSE)
  writeLines(c("",paste("CELLS",nrow(tets),nrow(tets)*5)),con=f)
  write.table(cbind(4,tets-1), file=f, row.names = FALSE, col.names=FALSE)
  writeLines(c("",paste("CELL_TYPES",nrow(tets))),con=f)
  write.table(cbind(rep(10,nrow(tets))), file=f, row.names = FALSE, col.names=FALSE)
  close(f)
  invisible(filename)
}