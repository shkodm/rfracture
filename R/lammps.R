#' Read and write balls from and to lammps data files
#' 
#' @param filename file path
#' 
#' @import utils
#' @export
read.lammps.data = function(filename) {
  f = file(filename,"r")
  l = readLines(f,11)
  l = strsplit(l," ")
  n = as.integer(l[[3]][1])
  tab = read.table(f,header=FALSE,nrows=n)
  l2 = readLines(f,3)
  if (l2[2] == "Velocities") {
    write.velocities = TRUE
    tabv = read.table(f,header=FALSE,nrows=n)
  } else {
    write.velocities = FALSE
    tabv = matrix(0, n, 7)
  }	
  close(f)
  B = data.frame(
    type=tab[,2],
    r=tab[,3]/2,
    x=tab[,5],
    y=tab[,6],
    z=tab[,7],
    dens = tab[,4],
    ux = tabv[,2],
    uy = tabv[,3],
    uz = tabv[,4],
    wx = tabv[,5],
    wy = tabv[,6],
    wz = tabv[,7]
  )
  lim = sapply(l[6:8],function(x) as.numeric(x[1:2]))
  list(B=B, xlim=lim[,1], ylim=lim[,2], zlim=lim[,3],write.velocities = write.velocities)
}

#' @rdname read.lammps.data
#' 
#' @param obj balls (list from read.lammps.data or data.frame)
#' @param filename file path
#' @param density density of balls exported
#' @param comment comment placed in the header
#' @param xlim,ylim,zlim extent of the domain exported in the lammps file
#' 
#' @import utils
#' @export
write.lammps.data = function(obj=list(B), B, file=stdout(), write.velocities, density=2.5, comment="R generated balls", xlim=range(obj$B$x), ylim=range(obj$B$y), zlim=range(obj$B$z)) {
  # atom-ID atom-type diameter density x y z
  if (! is.list(obj)) stop("obj in write.lammps.data has to be a list")
  for (n in c("xlim","ylim","zlim")) {
    mis = eval(substitute(missing(X), list(X = n)))
    if ((!mis) || (!n %in% names(obj))) obj[[n]] = get(n)
  }
  if ((!missing(density)) || (!"density" %in% names(obj))) obj$B$dens = density
  if (missing(write.velocities)) {
    if ("write.velocities" %in% names(obj)) {
      write.velocities = obj$write.velocities
    } else {
      write.velocities = any(c("ux","uy","uz") %in% names(obj$B))  
    }
  } 
  if ((!"density" %in% names(obj))) obj$B$dens = density
  if (is.character(file)) {
    f = file(file,"w")
  } else {
    f = file
  }
  writeLines(c(
    "LAMMPS data file via R",
    "",
    paste(nrow(obj$B),"atoms"),
    paste(max(obj$B$type), "atom","types"),
    "",
    paste(obj$xlim[1],obj$xlim[2],"xlo","xhi"),
    paste(obj$ylim[1],obj$ylim[2],"ylo","yhi"),
    paste(obj$zlim[1],obj$zlim[2],"zlo","zhi"),
    "",
    "Atoms",
    ""),f)
  tab = data.frame(
    seq_len(nrow(obj$B)),
    obj$B$type,
    obj$B$r*2,
    obj$B$dens,
    obj$B$x, 
    obj$B$y, 
    obj$B$z)
  write.table(tab, file=f,row.names=FALSE,col.names=FALSE)
  if (write.velocities) {
    for (n in c("ux","uy","uz","wx","wy","wz")) if (! n %in% names(obj$B)) obj$B[[n]] = 0
    tabv = data.frame(
      seq_len(nrow(obj$B)),
      obj$B$ux, 
      obj$B$uy, 
      obj$B$uz,
      obj$B$wx, 
      obj$B$wy, 
      obj$B$wz)
    writeLines(c(
      "",
      "Velocities",
      ""),f)
    write.table(tabv, file=f,row.names=FALSE,col.names=FALSE)
  }
  if (is.character(file)) {
    close(f)
  }
  obj
}
