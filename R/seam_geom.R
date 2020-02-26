
fracture.touching = function(obj,touch="exclude") {
  P = obj$points
  i = obj$triangles
  sel = P$f1 == P$f2

  sel = sel[i]
  dim(sel) = dim(i)
  sel = rowSums(sel) != 3
  if (touch == "exclude") {
    
  } else if (touch == "include") {
    sel[] = TRUE
  } else if (touch == "only") {
    sel = ! sel
  } else stop("invalid value for touch")
  i = i[sel,]
  sel = rep(FALSE,nrow(P))
  sel[i[,1]] = TRUE
  sel[i[,2]] = TRUE
  sel[i[,3]] = TRUE
  
  ni = rep(0,nrow(P))
  ni[sel] = 1:sum(sel)
  i[] = ni[i]
  i = i[rowSums(i == 0) == 0,]
  P = P[sel,]
  obj$points=P
  obj$triangles=i
  return(obj)
  
}

fracture.volume = function(obj) {
  v = obj$points[,c("x","y")]
  v1 = v[obj$triangles[,2],] - v[obj$triangles[,1],]
  v2 = v[obj$triangles[,3],] - v[obj$triangles[,1],]
  a = abs(1/2*(v1[,1]*v2[,2] - v1[,2]*v2[,1]))
  h = 1/3*(obj$points$h[obj$triangles[,1]] + obj$points$h[obj$triangles[,2]] + obj$points$h[obj$triangles[,3]])
  sum(a*h)
}

