#' Makes a 3d plot of the fracture
#'
#' @export
#' @import rgl
fracture3d = function(obj,type=c("top","bottom"),top="top" %in% type,bottom="bottom" %in% type,middle="middle" %in% type,col=c(2,3,4), add=FALSE) {
  if (length(col) == 1) col = rep(col,3)
  iv = as.vector(t(obj$triangles))
  if (!add) {
    clear3d()
    plot3d(diag(3))
    clear3d()
  }
  if (top) triangles3d(obj$points$f1[iv],obj$points$x[iv],obj$points$y[iv],col=col[1])
  if (bottom) triangles3d(obj$points$f2[iv],obj$points$x[iv],obj$points$y[iv],col=col[2])
  if (middle) triangles3d((obj$points$f1+obj$points$f2)[iv]/2,obj$points$x[iv],obj$points$y[iv],col=col[3])
}

border3d = function(obj, f1, f2, add=FALSE, ...) {
  edges = rbind(obj$triangles[,1:2],obj$triangles[,2:3],obj$triangles[,c(1,3)])
  sel = edges[,1] > edges[,2]
  edges[sel,] = edges[sel,c(2,1)]
  head(edges)
  edges = edges[order(edges[,1],edges[,2]),]
  a = !duplicated(edges)
  head(edges[a,])
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
    plot3d(diag(3))
    clear3d()
  }
  triangles3d(rbind(t1,t2),...)
}

#' Plot fracture matrix
#' 
#' @param obj fracture matrix to plot
#' @param field field to plot ("f1" - upper, "f2" - lower)
#' @param col.palette color palette to use for plotting
#' @param ... other options for plot (scatter plot)
#' 
#' @export
plot.fracturematrix = function(obj, field="f1", col.palette=c("black","red","yellow"), pch=16, cex=1, asp=1, ...) {
  if (! field %in% names(obj)) stop(field, "is not in obj")
  col = as.vector(obj[[field]])
  col = (col-min(col))/(max(col)-min(col))
  col = colorRamp(col.palette)(col)
  col = rgb(col,max=255)
  plot(obj$points[,1], obj$points[,2], col=col, pch=pch, cex=cex, asp=asp, ...)
  abline(h=-5:5)
  abline(v=-5:5)
  arrows(0, 0, obj$mat[,1], obj$mat[,2], angle = 15)
}

