library(fields)
source("rtri.R")

seam.balls = function(obj,
                      K = 4000,
                      Rmax = 0.01,
                      Rmin = 0.01,
                      B = data.frame(),
                      margin = 0.003,
                      margin_opt = margin/3,
                      relax_iterations = 10,
                      max_add = 200,
                      overshot = 100,
                      seed,mean.neighbor = 5,
                      iterations = ceiling(1.5*K/max_add),
                      dist = function(k,Rmin,Rmax) runif(k,Rmin,Rmax),
                      delete=TRUE) {
  P = obj$points
  i = obj$triangles
  T = data.frame(i1=i[,1],i2=i[,2],i3=i[,3])
  h = P$h[i]
  dim(h) = dim(i)
  T$h = rowMeans(h)
  
  sel =        P$x[T$i1] > 0 & P$x[T$i2] > 0 & P$x[T$i3] > 0
  sel = sel & (P$y[T$i1] > 0 & P$y[T$i2] > 0 & P$y[T$i3] > 0)
  sel = sel & (P$x[T$i1] < 1 & P$x[T$i2] < 1 & P$x[T$i3] < 1)
  sel = sel & (P$y[T$i1] < 1 & P$y[T$i2] < 1 & P$y[T$i3] < 1)
  T = T[sel,]

  
  ndel = 100
  dlog = NULL
  if (!missing(seed)) set.seed(seed)
  for (iteration in seq_len(iterations)) {
    if (K+ndel > nrow(B)) {
      k = min(K + overshot - nrow(B), max_add)
      p = pmax(0, T$h - 2 * (Rmin+margin))
      
      nB = data.frame(
        tri = sample.int(nrow(T), k, prob = p, replace = TRUE),
        r = dist(k,Rmin,Rmax)
      )
      nB$h1 = P$h[T$i1[nB$tri]] - 2 * (nB$r+margin)
      nB$h2 = P$h[T$i2[nB$tri]] - 2 * (nB$r+margin)
      nB$h3 = P$h[T$i3[nB$tri]] - 2 * (nB$r+margin)
      sel = nB$h1+nB$h2+nB$h3 > 0
      nB = nB[sel,,drop=FALSE]
      ret = rtri(nrow(nB),nB$h1,nB$h2,nB$h3)
      nB$y = P$x[T$i1[nB$tri]]*ret$w1 + P$x[T$i2[nB$tri]]*ret$w2 + P$x[T$i3[nB$tri]]*ret$w3
      nB$z = P$y[T$i1[nB$tri]]*ret$w1 + P$y[T$i2[nB$tri]]*ret$w2 + P$y[T$i3[nB$tri]]*ret$w3
      nB$f1 = P$f1[T$i1[nB$tri]]*ret$w1 + P$f1[T$i2[nB$tri]]*ret$w2 + P$f1[T$i3[nB$tri]]*ret$w3
      nB$f2 = P$f2[T$i1[nB$tri]]*ret$w1 + P$f2[T$i2[nB$tri]]*ret$w2 + P$f2[T$i3[nB$tri]]*ret$w3
      nB$x = nB$f2 + nB$r + margin + ret$h
      if (! (all(nB$x > nB$f2 + nB$r + margin) && all(nB$x < nB$f1 - nB$r - margin))) stop("balls do not fit in triangles")
      B = rbind(B, nB)
      cat("Created", nrow(nB), "balls\n")
    }
  
    n = nrow(B)
    Xi = c(1:n,1:n,1:n,1:n,rep(0,2*nrow(P)))
    N = n*4
    start_ds = Inf
    for (i in seq_len(relax_iterations + 1)) {
      X = cbind(
        c(B$x, B$x,   B$x, B$x,   P$f1, P$f2),
        c(B$y, B$y+1, B$y, B$y+1, P$x , P$x ),
        c(B$z, B$z, B$z+1, B$z+1, P$y , P$y ))
      #print(dim(X))
      ds = fields.rdist.near(X[1:N,], X[1:n,], delta = 2*(Rmax+margin_opt),mean.neighbor = mean.neighbor,max.points = mean.neighbor*n)
      tds = data.frame(d = ds$ra, i = ds$ind[,1], j = ds$ind[,2])
      ds = fields.rdist.near(X[(N+1):nrow(X),], X[1:n,], delta = Rmax+margin_opt, mean.neighbor = mean.neighbor, max.points = mean.neighbor*n)
      if (length(ds$ra) == 1) tds = rbind(tds, data.frame(d = ds$ra, i = ds$ind[1]+N, j = ds$ind[2]))
        else if (length(ds$ra) != 0) tds = rbind(tds, data.frame(d = ds$ra, i = ds$ind[,1]+N, j = ds$ind[,2]))
      sel = tds$i > tds$j
      tds = tds[sel,,drop=FALSE]
      tds$iball = Xi[tds$i]
      tds$v = tds$d - ifelse(tds$iball > 0, B$r[tds$iball], 0) - B$r[tds$j]
      #print(range(tds$v))
      if (i == 1) start_ds = { if (nrow(tds) > 0) min(tds$v) else Inf }
      if (all(tds$v > 0) || i == relax_iterations + 1) {
        finish_ds = { if (nrow(tds) > 0) min(tds$v) else Inf }
        cat("Changed the minimal distance from", start_ds, "to", finish_ds, "in", i-1, "iterations\n");
        break;
      }
      p = X[tds$i,,drop=FALSE] - X[tds$j,,drop=FALSE]
      pl = sqrt(rowSums(p^2))
      sel = pl < 1e-10
      if (any(sel)) {
        cat("Some zero distances\n")
        p = p[!sel,]; pl = pl[!sel]; tds = tds[!sel,]
      }
      p = p / pl * pmax(0,margin - tds$v)
      p = rbind(-p,p)
      pi = c(tds$j,tds$iball)
      sel = pi != 0
      p = p[sel,,drop=FALSE]
      pi = pi[sel]
      p = sapply(1:3, function(k) tapply(p[,k],pi,mean))
      pi = tapply(pi,pi,function(x)x[1])
      X[pi,] = X[pi,] + p
      #print(range(p))
      #print(range(pi))
      B$x = X[1:n,1]
      B$y = X[1:n,2] %% 1
      B$z = X[1:n,3] %% 1
    }
    #return (list(B,tds,X))
    if (delete) {
      delsel = NULL
      o = order(tds$v)
      tds = tds[o,]
      tds = tds[tds$v < 0,]
      while(nrow(tds) > 0) {
        sel = tds$j[1]
        delsel = c(delsel, sel)
        tds = tds[tds$iball != sel & tds$j != sel, , drop = FALSE]
      }
      cat("Deleted", length(delsel), "balls\n")
      if (length(delsel) > 0)  B = B[-delsel,]
    }
    dlog = rbind(dlog,data.frame(n=nrow(B)))
    plot(dlog$n)
    if (nrow(B) >= K) {
      cat("Finished.\n")
      break;
    }
  }
  B
}

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
