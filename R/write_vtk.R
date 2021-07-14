#' Writes a tet-mesh of the fracture geometry in a legacy VTK format
#' 
#' @param obj the fracture_geom object
#' @param filename name of the VTK file to write
#' @param ... arguments passed to extract.tet.mesh, such as transform, or genie.h
#' 
#' @export
write.vtk.tet.mesh = function(obj, filename, ...){
  ex = extract.tet.mesh(obj, type="interior", ...)
  f = file(filename,"w")
  writeLines(c(
    "# vtk DataFile Version 2.0",
    "Unstructured Grid",
    "ASCII",
    "DATASET UNSTRUCTURED_GRID",
    paste("POINTS",nrow(ex$points),"double")),con=f)
  write.table(ex$points, file=f, row.names = FALSE, col.names=FALSE)
  writeLines(c("",paste("CELLS",nrow(ex$tetrahedra),nrow(ex$tetrahedra)*5)),con=f)
  write.table(cbind(4,ex$tetrahedra[,1:4]-1), file=f, row.names = FALSE, col.names=FALSE)
  writeLines(c("",paste("CELL_TYPES",nrow(ex$tetrahedra))),con=f)
  write.table(cbind(rep(10,nrow(ex$tetrahedra))), file=f, row.names = FALSE, col.names=FALSE)
  close(f)
  invisible(filename)
}
