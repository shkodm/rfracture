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
#' ret = slice(ret, eps=1e-2, flatten="above")
#' clear3d()
#' wire3d(as.mesh3d(ret))
#' 
#' @export
slice = function(obj, by="h", flatten="edge", value="all", level=0, eps=1e-10) {
  # select points below level
  psel = obj$points[,by] < level
  
  #select triangles crossing level
  i = cbind(psel[obj$triangles[,1]],psel[obj$triangles[,2]],psel[obj$triangles[,3]])
  j = rowSums(i)
  tsel = j > 0 & j < 3
  tr = obj$triangles[tsel,]
  
  below = j == 3
  
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
  
  if (max(abs(obj$points[p_on_edge,by] - level)) > eps) stop("error in slice above ", eps)
  
  flat = rep(FALSE,nrow(obj$points))
  if (flatten == "edge") {
    flat[p_on_edge] = TRUE
  } else if (flatten == "above") {
    flat[p_on_edge] = TRUE
    flat = flat | obj$points[,by] >= level
  } else if (flatten == "below") {
    flat[p_on_edge] = TRUE
    flat = flat | obj$points[,by] <= level
  } else if (flatten == "none") {
  } else stop("Unknown flatten parameter")
  obj$points[flat,by] = level
  if (by == "h") {
    obj$points[flat,"f1"] = obj$points[flat,"fm"] + level/2
    obj$points[flat,"f2"] = obj$points[flat,"fm"] - level/2
  }

  id = matrix(id,ncol=2)

  ntr = NULL
  sel = !is.na(id[,1]) & !is.na(id[,2])
  cat(sum(sel),"triangles trimmed\n")
  ntr = rbind(ntr,
    cbind(tr[sel,1],id[sel,1],id[sel,2]),
    cbind(tr[sel,2],id[sel,2],id[sel,1]),
    cbind(tr[sel,2],tr[sel,3],id[sel,2]))
  sel = is.na(id[,1]) & !is.na(id[,2])
  cat(sum(sel),"triangles cut through 1st vertex\n")
  ntr = rbind(ntr,
              cbind(tr[sel,1],tr[sel,2],id[sel,2]),
              cbind(tr[sel,2],tr[sel,3],id[sel,2]))
  sel = !is.na(id[,1]) & is.na(id[,2])
  cat(sum(sel),"triangles cut through 2nd vertex\n")
  ntr = rbind(ntr,
              cbind(tr[sel,1],id[sel,1],tr[sel,3]),
              cbind(tr[sel,2],tr[sel,3],id[sel,1]))
  sel = is.na(id[,1]) & is.na(id[,2])
  cat(sum(sel),"triangles cut on edge\n")
  ntr = rbind(ntr,
              cbind(tr[sel,1],tr[sel,2],tr[sel,3]))


  obj$triangles = rbind(obj$triangles[!tsel,], ntr)

  return(obj)
}
