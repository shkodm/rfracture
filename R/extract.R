sortRows = function(M) {
  if (! is.matrix(M)) stop("'M' is not a matrix")
  matrix(M[order(row(M), M)], ncol=ncol(M), byrow = TRUE)
}

#' Extracts points and triangles from fracture_geom
#' 
#' @param obj the fracture_geom object
#' @param type selects one or more of: top, bottom, middle, border, edge
#' @param tranform the function used for transform the points
#' 
#' @export
extract.tet.mesh = function(obj, type=c("top","bottom"), transform=function(points) points, genie.h) {
  if (! "fracture_geom" %in% class(obj)) stop("'obj' not of fracture_geom class")
  if (! missing(type)) if (! is.character(type)) stop("'type' not a string")
  points = data.frame(
    x = c(obj$points[, "x"], obj$points[, "x"], obj$points[, "x"]),
    y = c(obj$points[, "y"], obj$points[, "y"], obj$points[, "y"]),
    z = c(obj$points[,"f1"], obj$points[,"f2"], obj$points[,"fm"])
  )
  k = nrow(obj$points)
  simplexes = NULL

  add.simplexes = function(v1,v2=NA, v3=NA, v4=NA, tag) {
    simplexes <<- rbind(simplexes, data.frame(v1=v1,v2=v2,v3=v3,v4=v4,tag=tag,stringsAsFactors=TRUE))
  }
  
  ed = obj$edge
  tr = obj$triangles
  br = ifelse(obj$border, "border","edge")
  
  add.simplexes( v1 = ed[,1],     v2 = ed[,2],      tag = paste("top",br))
  add.simplexes( v1 = ed[,1]+k,   v2 = ed[,2]+k,    tag = paste("bottom",br))
  add.simplexes( v1 = ed[,1]+2*k, v2 = ed[,2]+2*k,  tag = paste("middle",br))
  ed = sortRows(ed)
  add.simplexes( v1 = tr[,1],     v2 = tr[,2],      v3 = tr[,3],      tag = "top")
  add.simplexes( v1 = tr[,1]+k,   v2 = tr[,2]+k,    v3 = tr[,3]+k,    tag = "bottom")
  add.simplexes( v1 = tr[,1]+2*k, v2 = tr[,2]+2*k,  v3 = tr[,3]+2*k,  tag = "middle")
  add.simplexes( v1 = ed[,2],     v2 = ed[,2]+k,    v3 = ed[,1],      tag = br)
  add.simplexes( v1 = ed[,1]+k,   v2 = ed[,1],      v3 = ed[,2]+k,    tag = br)
  tr = sortRows(tr)
  add.simplexes( v1 = tr[,1]+k,   v2 = tr[,1],      v3 = tr[,2],      v4 = tr[,3],     tag = "interior")
  add.simplexes( v1 = tr[,2]+k,   v2 = tr[,2],      v3 = tr[,3],      v4 = tr[,1]+k,   tag = "interior")
  add.simplexes( v1 = tr[,3]+k,   v2 = tr[,3],      v3 = tr[,1]+k,    v4 = tr[,2]+k,   tag = "interior")
  
  simplexes$order = ifelse(is.na(simplexes$v4),ifelse(is.na(simplexes$v3),ifelse(is.na(simplexes$v2),ifelse(is.na(simplexes$v1),0,1),2),3),4)
  if (! missing(type)) simplexes = simplexes[simplexes$tag %in% type, , drop=FALSE]
  selp = na.omit(unique(c(simplexes$v1,simplexes$v2,simplexes$v3,simplexes$v4)))
  map = rep(0,nrow(points)); map[selp] = seq_len(length(selp))
  
  if (!missing(genie.h)) {
    cl = genieclust::gclust(points)
    map = cutree(cl, h = genie.h)
    points = aggregate(points, by=list(map), mean)
    points$Group.1 = NULL
    simplexes$v1 = map[simplexes$v1]
    simplexes$v2 = map[simplexes$v2]
    simplexes$v3 = map[simplexes$v3]
    simplexes$v4 = map[simplexes$v4]
  }
  
  dups = pmax(simplexes$v1 == simplexes$v2,simplexes$v1 == simplexes$v3,simplexes$v2 == simplexes$v3,simplexes$v1 == simplexes$v4,simplexes$v2 == simplexes$v4,simplexes$v3 == simplexes$v4,na.rm = TRUE)
  simplexes = simplexes[dups == 0,]
  
  points = points[selp,,drop=FALSE]
  simplexes$v1 = map[simplexes$v1]
  simplexes$v2 = map[simplexes$v2]
  simplexes$v3 = map[simplexes$v3]
  simplexes$v4 = map[simplexes$v4]
  
  points = transform(points)
  list(
    points     = points,
    vertexes   = simplexes[simplexes$order == 1,,drop=FALSE],
    edges      = simplexes[simplexes$order == 2,,drop=FALSE],
    triangles  = simplexes[simplexes$order == 3,,drop=FALSE],
    tetrahedra = simplexes[simplexes$order == 4,,drop=FALSE]
  )
}
