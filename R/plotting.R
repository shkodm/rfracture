#' Makes a 3d plot of the fracture
#' 
#' @param obj the fracture_geom object to plot in 3d
#' @param type the surface to plot: top, bottom, middle
#' @param col vector of 3 colors used for top, bottom and middle surface
#' @param add if TRUE, add plot to existing rgl window
#' 
#' @import rgl
#' @export
fracture3d = function(obj, type=c("top","bottom"), col=c(2,3,4), add=FALSE) {
  if (length(col) == 1) col = rep(col,3)
  iv = as.vector(t(obj$triangles))
  if (!add) {
    clear3d()
    aspect3d("iso")
  }
  if ("top"    %in% type) triangles3d(obj$points$f1[iv],obj$points$x[iv],obj$points$y[iv],col=col[1])
  if ("bottom" %in% type) triangles3d(obj$points$f2[iv],obj$points$x[iv],obj$points$y[iv],col=col[2])
  if ("middle" %in% type) triangles3d((obj$points$f1+obj$points$f2)[iv]/2,obj$points$x[iv],obj$points$y[iv],col=col[3])
}

#' Makes a mesh3d object from fracture
#' 
#' @param obj the fracture_geom object to plot in 3d
#' @param type the surface to plot: top, bottom, middle
#' @param ... other arguments passed to tmesh3d
#' 
#' @import rgl
#' @export
as.mesh3d.fracture_geom = function(obj, type="top", ...) {
  if (type == "top") {
    P = obj$points[,c("f1","x","y")]
  } else if (type == "bottom") {
    P = obj$points[,c("f2","x","y")]
  } else if (type == "middle") {
    P = obj$points[,c("fm","x","y")]
  } else stop("unknown type:", type)
  tmesh3d(vertices = t(as.matrix(P)), indices = t(as.matrix(obj$triangles)),homogeneous = FALSE, ...)
}


border3d = function(obj, f1, f2, add=FALSE, ...) {
  edges = rbind(obj$triangles[,1:2],obj$triangles[,2:3],obj$triangles[,c(1,3)])
  sel = edges[,1] > edges[,2]
  edges[sel,] = edges[sel,c(2,1)]
  edges = edges[order(edges[,1],edges[,2]),]
  a = !duplicated(edges)
  sel = which(table(cumsum(a)) == 1)
  edges = edges[a,][sel,]
  
  plot(obj$points$x[edges[,1]],obj$points$y[edges[,2]],asp=1)
  
  if (missing(f1)) f1 = obj$points$f1
  if (missing(f2)) f2 = obj$points$f2
  if (length(f1) == 1) f1 = rep(f1, nrow(obj$points))
  if (length(f2) == 1) f2 = rep(f2, nrow(obj$points))
  open1 = f1[edges[,1]] != f2[edges[,1]]
  open2 = f1[edges[,2]] != f2[edges[,2]]
  
  t1 = rbind(
    f1[edges[open1,1]],obj$points$x[edges[open1,1]],obj$points$y[edges[open1,1]],
    f2[edges[open1,1]],obj$points$x[edges[open1,1]],obj$points$y[edges[open1,1]],
    f2[edges[open1,2]],obj$points$x[edges[open1,2]],obj$points$y[edges[open1,2]])
  dim(t1) = c(3,sum(open1)*3)
  t1 = t(t1)
  t2 = rbind(
    f2[edges[open2,2]],obj$points$x[edges[open2,2]],obj$points$y[edges[open2,2]],
    f1[edges[open2,2]],obj$points$x[edges[open2,2]],obj$points$y[edges[open2,2]],
    f1[edges[open2,1]],obj$points$x[edges[open2,1]],obj$points$y[edges[open2,1]])
  dim(t2) = c(3,sum(open2)*3)
  t2 = t(t2)
  if (!add) {
    clear3d()
    aspect3d("iso")
  }
  triangles3d(rbind(t1,t2),...)
}

#' Plot fracture matrix
#' 
#' @param x fracture matrix to plot
#' @param field field to plot ("f1" - upper, "f2" - lower)
#' @param col.palette color palette to use for plotting
#' @param pch,cex,asp,... other options for plot (scatter plot)
#' 
#' @import graphics
#' @import grDevices
#' @export
plot.fracture_matrix = function(x, field="f1", col.palette=c("black","red","yellow"), pch=16, cex=1, asp=1, ...) {
  obj = x
  if (! field %in% names(obj)) stop(field, "is not in obj")
  if (length(obj$dims) == 1) {
    matplot(obj$points,cbind(obj$f1,obj$f2), lty=1, type="l", asp=asp, ...)
  } else if (length(obj$dims) == 2) {
    col = as.vector(obj[[field]])
    col = (col-min(col))/(max(col)-min(col))
    col = colorRamp(col.palette)(col)
    col = rgb(col,maxColorValue=255)
    plot(obj$points[,1], obj$points[,2], col=col, pch=pch, cex=cex, asp=asp, ...)
    arrows(0, 0, obj$span[,1], obj$span[,2], angle = 15)
    seg = function(a,b) segments(a[,1],a[,2],b[,1],b[,2])
    seg(cbind(-1,0:5) %*% obj$period, cbind(6,0:5) %*% obj$period)
    seg(cbind(0:5,-1) %*% obj$period, cbind(0:5,6) %*% obj$period)
  } else stop("Plot not implemented with highdimensional fractures")
}

