#' Slices the fracture geometry
#' 
#' Slices the fracture with respect to a specific property, and threshold level
#' 
#' @param obj fracture_geom object
#' @param by name of the property/field by which to slice
#' @param level threshold level on the property
#' @param value type of the returned object (see return value)
#' @param flatten if to snap values to level
#' @param eps the numerical accuracy
#' 
#' @return
#' If value is "all", the function returns the same fracture_geom object, with some of the triangles sliced.
#' If value is "above", the function returns only the part of the geometry, for which property > level
#' "below" returns the opposite.
#' "equal" returns a list(points, edges), where points are points with property == level, and edges are the connecting edges (two columns)
#' 
#' @examples
#' library(rgl)
#' ret = fracture_geom(refine=4)
#' ret2 = slice(ret, eps=1e-2, value="above")
#' clear3d()
#' wire3d(as.mesh3d(ret2),col=2)
#' ret2 = slice(ret, eps=1e-2, value="below")
#' wire3d(as.mesh3d(ret2),col=3)
#' 
#' @export
slice = function(obj, by="h", flatten="edge", value="all", level=0, eps=1e-10, verbose=FALSE) {
  # select points below level
  psel = obj$points[,by] < level
  
  #select triangles crossing level
  i = cbind(psel[obj$triangles[,1]],psel[obj$triangles[,2]],psel[obj$triangles[,3]])
  j = rowSums(i)
  tsel = j > 0 & j < 3
  tr = obj$triangles[tsel,]
  
  #sort vertices in triangles so that first is on the other side then other two
  i = i[tsel,]
  j = j[tsel] == 1
  tr = rbind(tr[i[,1] == j,c(1,2,3)],tr[i[,2] == j,c(2,3,1)],tr[i[,3] == j,c(3,1,2)])

  #generate unique indexes
  i = cbind(
    c(tr[,1],tr[,1]),
    c(tr[,2],tr[,3]))
  i[i[,1] > i[,2], ] = i[i[,1] > i[,2], 2:1]
  id = paste(i[,1],i[,2],sep="_")

  a = obj$points[i[,1],by] - level
  b = obj$points[i[,2],by] - level
  t = a/(a-b)
  t[a-b == 0] = 0.5
  t[t<eps] = 0
  t[t>1-eps] = 1
  
  np = (1-t) * obj$points[i[,1],] + obj$points[i[,2],]*t
  id[t == 0] = NA
  id[t == 1] = NA
  id = as.integer(factor(id)) + nrow(obj$points)
  
  npsel = ! (duplicated(id) | is.na(id)) 
  np = np[npsel,]
  np_id = id[npsel]
  np = np[order(np_id),]
  np_id = np_id[order(np_id)]
  
  p_on_edge = c(i[t == 0,1], i[t == 1, 2], np_id)
  obj$points = rbind(obj$points,np)
  
  on_edge = rep(FALSE,nrow(obj$points))
  on_edge[p_on_edge] = TRUE
  above = obj$points[,by] >= level & !on_edge
  below = obj$points[,by] <= level & !on_edge
  
  if (max(abs(obj$points[p_on_edge,by] - level)) > eps) stop("error in slice above ", eps)
  
  
  if (flatten == "edge") {
    flat = on_edge
  } else if (flatten == "above") {
    flat = on_edge | above
  } else if (flatten == "below") {
    flat = on_edge | below
  } else if (flatten == "none") {
    flat = rep(FALSE,nrow(obj$points))
  } else stop("Unknown flatten parameter")
  obj$points[flat,by] = level
  if (by == "h") {
    obj$points[flat,"f1"] = obj$points[flat,"fm"] + level/2
    obj$points[flat,"f2"] = obj$points[flat,"fm"] - level/2
  }
  id = matrix(id,ncol=2)

  ntr = NULL
  edge = NULL
  sel = !is.na(id[,1]) & !is.na(id[,2])
  if (verbose) cat(sum(sel),"triangles trimmed\n")
  ntr = rbind(ntr,
    cbind(tr[sel,1],id[sel,1],id[sel,2]),
    cbind(tr[sel,2],id[sel,2],id[sel,1]),
    cbind(tr[sel,2],tr[sel,3],id[sel,2]))
  sel = is.na(id[,1]) & !is.na(id[,2])
  if (verbose) cat(sum(sel),"triangles cut through 1st vertex\n")
  ntr = rbind(ntr,
              cbind(tr[sel,1],tr[sel,2],id[sel,2]),
              cbind(tr[sel,2],tr[sel,3],id[sel,2]))
  sel = !is.na(id[,1]) & is.na(id[,2])
  if (verbose) cat(sum(sel),"triangles cut through 2nd vertex\n")
  ntr = rbind(ntr,
              cbind(tr[sel,1],id[sel,1],tr[sel,3]),
              cbind(tr[sel,2],tr[sel,3],id[sel,1]))
  sel = is.na(id[,1]) & is.na(id[,2])
  if (verbose) cat(sum(sel),"triangles cut on edge\n")
  ntr = rbind(ntr,
              cbind(tr[sel,1],tr[sel,2],tr[sel,3]))

  sel = (id[,1] != id[,2])
  edge = rbind(edge, cbind(id[sel,1],id[sel,2]))

  obj$triangles = rbind(obj$triangles[!tsel,], ntr)

  if (value == "edge") {
    psel = on_edge
  } else if (value == "above") {
    psel = on_edge | above
  } else if (value == "below") {
    psel = on_edge | below
  } else if (value == "all") {
    psel = rep(TRUE,nrow(obj$points))
  } else if (value == "none") {
    psel = rep(FALSE,nrow(obj$points))
  } else stop("Unknown flatten parameter")
  
  
  tsel = psel[obj$triangles]
  dim(tsel) = dim(obj$triangles)
  tsel = rowSums(tsel) == 3
  obj$triangles = obj$triangles[tsel,]
  esel = psel[edge]
  dim(esel) = dim(edge)
  
  psel = rep(FALSE,nrow(obj$points))
  psel[obj$triangles[,1]] = TRUE
  psel[obj$triangles[,2]] = TRUE
  psel[obj$triangles[,3]] = TRUE
  psel[edge[,1]] = TRUE
  psel[edge[,2]] = TRUE
  
  ni = rep(0,nrow(obj$points))
  ni[psel] = 1:sum(psel)
  obj$triangles[] = ni[obj$triangles]
  edge[] = ni[edge]
  if (any(rowSums(obj$triangles == 0) != 0)) stop("what?")
  obj$points = obj$points[psel,]
  obj$edge = edge

  return(obj)
}
