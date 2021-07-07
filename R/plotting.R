#' Makes a 3d plot of the fracture
#' 
#' @param obj the fracture_geom object to plot in 3d
#' @param type the surface to plot: top, bottom, middle
#' @param col vector of 3 colors used for top, bottom and middle surface
#' @param add if TRUE, add plot to existing rgl window
#' 
#' @export
fracture3d = function(obj, type=c("top","bottom"), col, edge.col=NA, vertex.col=NA, asp="iso", add=FALSE, ...) {
  ex = extract.tet.mesh(obj, type=type)
  n = nlevels(ex$triangles$tag)
  if (missing(col)) col = seq_len(n)+1
  if (length(col) == 1) col = rep(col, n)
  if (length(col) != n) stop("need colors for: ",paste(levels(ex$triangles$tag),collapse=", "))
  n = nlevels(ex$edges$tag)
  if (missing(edge.col)) edge.col = seq_len(n)+1
  if (length(edge.col) == 1) edge.col = rep(edge.col, n)
  if (length(edge.col) != n) stop("need colors for: ",paste(levels(ex$edges$tag),collapse=", "))
  
  if (!add) {
    rgl::plot3d(ex$points, type="n", asp=asp, ...)
  }
  ie = as.vector(t(ex$edges[,1:2]))
  it = as.vector(t(ex$triangles[,1:3]))
  if (length(ie) != 0) rgl::segments3d(ex$points[ie,c("x","y","z")], col=rep(edge.col[ex$edges$tag],each=2))
  if (length(it) != 0) rgl::triangles3d(ex$points[it,c("x","y","z")],col=rep(col[ex$triangles$tag],each=3))
  NULL
}

#' Makes a mesh3d object from fracture
#' 
#' @param obj the fracture_geom object to plot in 3d
#' @param type the surface to plot: top, bottom, middle
#' @param ... other arguments passed to tmesh3d
#' 
#' @export
as.mesh3d.fracture_geom = function(obj, type=c("middle","top","bottom",""), ...) {
  if (type == "top") {
    P = obj$points[,c("f1","x","y")]
  } else if (type == "bottom") {
    P = obj$points[,c("f2","x","y")]
  } else if (type == "middle") {
    P = obj$points[,c("fm","x","y")]
  } else stop("unknown type:", type)
  obj = list(
    vb = rbind(t(as.matrix(P)),1),
    it = t(as.matrix(obj$triangles)),
    material = list(),
    normals = NULL,
    meshColor = "vertices",
    texcoords = NULL
  )
  class(obj) = c("mesh3d", "shape3d")
  obj
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

#' Writes an stl file
#' 
#' @export
write.stl = function (x, ...) UseMethod("write.stl")

#' Writes a mesh3d object (or a list) to an stl file
#' 
#' @param x list of meshes to write
#' @param con file connection of filename to write to
#' @param ascii if TRUE, write in ASCII format (discouraged)
#' 
#' @rdname write.stl
#' @export
write.stl.default = function(x,con,ascii=FALSE) {
  if (all(class(x) == "list")) {
    meshes = x
  } else {
    meshes = list(x)
  }
  meshes = lapply(meshes,function(mesh) {
    if ("mesh3d" %in% class(mesh)) {
      mesh
    } else if (requireNamespace("rgl", quietly = TRUE)) {
      rgl::as.mesh3d(mesh)
    } else {
      NULL
    }
  })
  asEuclidean = function (x) {
    if (is.matrix(x) && dim(x)[2] == 4) 
      return(x[, 1:3]/x[, 4])
    else if (length(x) == 4) 
      return(c(x[1]/x[4], x[2]/x[4], x[3]/x[4]))
    else stop("'x' is not row vectors(s)")
  }
  writeHeader <- function() {
    ident <- paste(filename, " produced by RGL\n")
    if (ascii) 
      cat("solid ", ident, file = con)
    else {
      padding <- paste(rep(" ", 80), collapse = "")
      ident <- substr(paste("binary", ident, padding), 
                      1, 80)
      writeChar(ident, con, nchars = 80, useBytes = TRUE, 
                eos = NULL)
      writeBin(0L, con, size = 4, endian = "little")
    }
  }
  xprod = function(v1,v2) cbind(v1[2]*v2[3]-v2[2]*v1[3], v1[3]*v2[1]-v2[3]*v1[1], v1[1]*v2[2]-v2[1]*v1[2])
  normalize = function(v) v/sqrt(sum(v^2))
  triangles <- 0
  writeTriangles <- function(vertices) {
    if (nrow(vertices)%%3 != 0) 
      stop("Need 3N vertices")
    n <- nrow(vertices)/3
    for (i in seq_len(n)) {
      vec0 <- vertices[3 * i - 2, ]
      vec1 <- vertices[3 * i - 1, ]
      vec2 <- vertices[3 * i, ]
      normal <- normalize(xprod(vec2 - vec0, vec1 - vec0))
      if (ascii) {
        cat("facet normal ", normal, "\n", file = con)
        cat("outer loop\n", file = con)
        cat("vertex ", vec0, "\n", file = con)
        cat("vertex ", vec1, "\n", file = con)
        cat("vertex ", vec2, "\n", file = con)
        cat("endloop\n", file = con)
        cat("endfacet\n", file = con)
      }
      else {
        writeBin(c(normal, vec0, vec1, vec2), con, size = 4, 
                 endian = "little")
        writeBin(0L, con, size = 2, endian = "little")
      }
    }
    triangles <<- triangles + n
  }
  writeQuads <- function(vertices) {
    if (nrow(vertices)%%4 != 0) 
      stop("Need 4N vertices")
    n <- nrow(vertices)/4
    indices <- rep(seq_len(n) * 4, each = 6) - rep(c(3, 2, 
                                                     1, 3, 1, 0), n)
    writeTriangles(vertices[indices, ])
  }
  writeSurface <- function(id) {
    vertices <- rgl.attrib(id, "vertices")
    dims <- rgl.attrib(id, "dim")
    nx <- dims[1]
    nz <- dims[2]
    indices <- integer(0)
    for (j in seq_len(nx - 1) - 1) {
      v1 <- j + nx * (seq_len(nz) - 1) + 1
      indices <- c(indices, rbind(v1[-nz], v1[-nz] + 1, 
                                  v1[-1] + 1, v1[-1]))
    }
    writeQuads(vertices[indices, ])
  }
  writeMesh <- function(mesh, scale = 1, offset = c(0, 0, 0)) {
    vertices <- asEuclidean(t(mesh$vb)) * scale
    vertices <- vertices + rep(offset, each = nrow(vertices))
    nt <- length(mesh$it)/3
    nq <- length(mesh$ib)/4
    newverts <- matrix(NA, 3 * nt + 6 * nq, 3)
    if (nt) 
      newverts[1:(3 * nt), ] <- vertices[c(mesh$it), ]
    if (nq) 
      newverts[3 * nt + 1:(6 * nq), ] <- vertices[c(mesh$ib[1:3, 
                                                            ], mesh$ib[c(1, 3, 4), ]), ]
    writeTriangles(newverts)
  }
  avgScale <- function() {
    bbox <- par3d("bbox")
    ranges <- c(bbox[2] - bbox[1], bbox[4] - bbox[3], bbox[6] - 
                  bbox[5])
    if (prod(ranges) == 0) 
      1
    else exp(mean(log(ranges)))
  }
  if (is.character(con)) {
    con <- file(con, if (ascii) 
      "w"
      else "wb")
    on.exit(close(con))
  }
  filename <- summary(con)$description
  writeHeader()
  for (mesh in meshes) writeMesh(mesh)
  if (!ascii) {
    seek(con, 80)
    writeBin(as.integer(triangles), con, size = 4, endian = "little")
  }
  invisible(filename)
}

#' Writes a mesh3d object (or a list) to an stl file
#' 
#' @param x list of meshes to write
#' @param con file connection of filename to write to
#' @param ascii if TRUE, write in ASCII format (discouraged)
#' 
#' @export
write.stl.fracture_geom = function(x,con,ascii=FALSE, type=c("top","bottom"), ...) {
  meshes = lapply(type, function(type) as.mesh3d.fracture_geom(x, type, ...))
  write.stl(meshes, con=con, ascii=ascii)
}