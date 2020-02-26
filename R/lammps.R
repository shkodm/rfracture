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
  close(f)
  dens = tab[,4]
  B = data.frame(r=tab[,3]/2, x=tab[,5], y=tab[,6], z=tab[,7])
  lim = sapply(l[6:8],function(x) as.numeric(x[1:2]))
  list(B=B, dens, xlim=lim[,1], ylim=lim[,2], zlim=lim[,3])
}

#' @rdname read.lammps.data
#' 
#' @param B balls (data.frame)
#' @param filename file path
#' @param density density of balls exported
#' @param comment comment placed in the header
#' @param xlim,ylim,zlim extent of the domain exported in the lammps file
#' 
#' @import utils
#' @export
write.lammps.data = function(B, filename, density=2, comment="R generated balls", xlim=range(B$x), ylim=range(B$y), zlim=range(B$z)) {
  # atom-ID atom-type diameter density x y z
  atoms = data.frame(
    id = 1:nrow(B),
    type = 1,
    diam = B$r*2,
    dens = density,
    x = B$x,
    y = B$y,
    z = B$z
  )
  f = file(filename, open = "w")
  cat(comment, "\n",file=f)
  cat("\n",file=f)
  cat(nrow(atoms),"atoms\n",file=f)
  cat("1 atom types\n",file=f)
  cat("0 bond types\n",file=f)
  cat("0 angle types\n",file=f)
  cat("\n",file=f)
  cat(xlim[1],xlim[2],"xlo xhi\n",file=f)
  cat(ylim[1],ylim[2],"ylo yhi\n",file=f)
  cat(zlim[1],zlim[2],"zlo zhi\n",file=f)
  cat("\n",file=f)
  cat("Atoms\n",file=f)
  cat("\n",file=f)
  write.table(atoms,file=f,row.names=FALSE,col.names=FALSE)
  close(f)
}
