sortRows = function(M) {
  if (! is.matrix(M)) stop("'M' is not a matrix")
  matrix(M[order(row(M), M)], ncol=ncol(M), byrow = TRUE)
}

#' Extracts points and triangles from fracture_geom
#' 
#' @param obj the fracture_geom object
#' @param type selects one or more of: top, bottom, middle, border, edge
#' 
#' @export
extract.tet.mesh = function(obj, type=c("top","bottom")) {
  if (! "fracture_geom" %in% class(obj)) stop("'obj' not of fracture_geom class")
  if (! is.character(type)) stop("'type' not a string")
  points = data.frame(
    x = c(obj$points[, "x"], obj$points[, "x"], obj$points[, "x"]),
    y = c(obj$points[, "y"], obj$points[, "y"], obj$points[, "y"]),
    z = c(obj$points[,"f1"], obj$points[,"f2"], obj$points[,"fm"])
  )
  k = nrow(obj$points)
  edges = NULL
  triangles = NULL
  tetrahedra = NULL
  
  ed = obj$edge
  tr = obj$triangles
  br = ifelse(obj$border, "border","edge")
  edges = data.frame(
    v1 = c(ed[,1],ed[,1]+k,ed[,1]+2*k),
    v2 = c(ed[,2],ed[,2]+k,ed[,2]+2*k),
    tag = c(paste("top",br),paste("bottom",br),paste("middle",br)),
    stringsAsFactors = TRUE
  )
  
  ed = sortRows(ed)
  triangles = data.frame(
    v1 = c(tr[,1], tr[,1]+k, tr[,1]+2*k, ed[,1], ed[,2]+k),
    v2 = c(tr[,2], tr[,2]+k, tr[,2]+2*k, ed[,2], ed[,1]+k),
    v3 = c(tr[,3], tr[,3]+k, tr[,3]+2*k, ed[,2]+k, ed[,1]),
    tag = c(
      rep("top",nrow(tr)),
      rep("bottom",nrow(tr)),
      rep("middle",nrow(tr)),
      br,
      br
    ),
    stringsAsFactors = TRUE
  )
  tr = sortRows(tr)
  tetrahedra = data.frame(
    v1 = c(tr[,1]+k,tr[,2]+k,tr[,3]+k),
    v2 = c(tr[,1]  ,tr[,2]  ,tr[,3]  ),
    v3 = c(tr[,2]  ,tr[,3]  ,tr[,1]+k),
    v4 = c(tr[,3]  ,tr[,1]+k,tr[,2]+k),
    tag = "interior",
    stringsAsFactors = TRUE
  )

  edges = edges[edges$tag %in% type,,drop=FALSE]
  triangles = triangles[triangles$tag %in% type,,drop=FALSE]
  tetrahedra = tetrahedra[tetrahedra$tag %in% type,,drop=FALSE]
  
  selp = unique(c(
    edges$v1, edges$v2,
    triangles$v1, triangles$v2, triangles$v3,
    tetrahedra$v1, tetrahedra$v2, tetrahedra$v3, tetrahedra$v4
  ))
  map = rep(0,nrow(points)); map[selp] = seq_len(length(selp))
  points = points[selp,,drop=FALSE]
  edges$v1      = map[edges$v1]
  edges$v2      = map[edges$v2]
  triangles$v1  = map[triangles$v1]
  triangles$v2  = map[triangles$v2]
  triangles$v3  = map[triangles$v3]
  tetrahedra$v1 = map[tetrahedra$v1]
  tetrahedra$v2 = map[tetrahedra$v2]
  tetrahedra$v3 = map[tetrahedra$v3]
  tetrahedra$v4 = map[tetrahedra$v4]
  
  list(
    points=points,
    edges=edges,
    triangles=triangles,
    tetrahedra=tetrahedra
  )
}
