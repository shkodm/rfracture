ordered_rnorm_mat = function(N, M, k=1, seed) {
  if (! missing(seed)) set.seed(seed)
  RN = t(matrix(rnorm(N*M*k*2),k*2,N*M))
  R = matrix(0+0i,N,M)
  circ = function(x,n) ifelse(x > n/2,x-n,x)
  fpZ = expand.grid(i=circ(1:N-1,N)/N,j=circ(1:M-1,M)/M)
  o = order(pmax(abs(fpZ$i),abs(fpZ$j)),atan2(fpZ$i,fpZ$j))
  lapply(1:k-1, function(i) {
    R[o] = RN[,2*i+1] + 1i*RN[,2*i+2]
    R
  })
}


seam.geom.my = function(totalScale = 0.2,
                       surfaceScale = 0.13,
                       centerlineScale = 0.0,
                       centerlineRadius = totalScale/2,
                       centerlineMean = totalScale/2,
                       gapScale = 0,
                       gapRadius = 0.05,
                       gapMean = totalScale - surfaceScale,
                       surfaceRadius1 = 0.05,
                       surfaceRadius2 = 0.01,
                       surfaceRatio = 0.90,
                       shape1=2,
                       shape2=2) {
  function(refine,seed) {
    bonds = c(centerlineMean - centerlineScale/2 - (gapMean+gapScale/2)/2 - surfaceScale/2,
              centerlineMean + centerlineScale/2 + (gapMean+gapScale/2)/2 + surfaceScale/2)
    cat("Max   bonds:",bonds [1],bonds [2],"\n")
    n = 6 * refine
    m = 5 * refine
    
    N = 5*n
    M = m
    
    R = matrix(0,N,M)
    X = row(R)-1
    Y = col(R)-1
    B = Y/m + X/m / 2
    A = X/n
    A1 = A %% 1
    B1 = B %% 1
    
    D = sqrt((A1-0.5)^2 + (B1-0.5)^2)
    
    conv.exp = function(d, scale=1) {
      K = exp(-D/d)
      scale * K/sqrt(sum(K^2))
    }
    
    RN = ordered_rnorm_mat(N,M,3,seed=seed)
    
    beta_random_field = function(k, scale=1, radius=1, mean=0, skew=0, shape1=2, shape2=shape1, K) {
      if (missing (K)) K = conv.exp(radius) else K = K/sqrt(sum(K^2))
      f = fft(RN[[k]]*fft(K),inverse = TRUE)/sqrt(length(K))
      move = function(f) {
        f = qbeta(pnorm(f+skew), shape1 = shape1, shape2=shape2)
        (f-0.5)*scale+mean
      }
      list(Re=move(Re(f)),Im=move(Im(f)))
    }
    
    
    centerlineField = beta_random_field(1,
                                        scale = centerlineScale,
                                        radius = centerlineRadius,
                                        mean = centerlineMean,
                                        shape1=shape1, shape2=shape2)
    gapField = beta_random_field(2,
                                 scale = gapScale,
                                 radius = gapRadius,
                                 mean = gapMean,
                                 shape1=shape1, shape2=shape2)
    surfaceField = beta_random_field(3,
                                     scale = surfaceScale,
                                     K=conv.exp(surfaceRadius1,surfaceRatio)+conv.exp(surfaceRadius2,1-surfaceRatio),
                                     shape1=shape1, shape2=shape2)
    
    
    list(f1 = centerlineField[[1]] + gapField[[1]]/2 + surfaceField[[1]],
      f2 = centerlineField[[1]] - gapField[[1]]/2 + surfaceField[[2]])
  }
}

# profile from the article: Fluid flow through rough fractures in rocks. II: A new matchingmodel for rough rock fractures
ogilve.corr.profile = function(ML=0.5, TL=0, MinMF=0, MaxMF=1) function(lambda) {
  ifelse(lambda >= ML+TL/2, 1, ifelse(lambda <= ML-TL/2, 0, (ML-TL/2-lambda)*(lambda-(ML+3*TL/2))/(TL*TL) ))*(MaxMF-MinMF)+MinMF
}

# spectral density based on Brown: Simple mathematical model of a rough fracture
seam.geom.brown = function(scale=1, alpha=2, corr.profile, closed=0.1, gap=NA)  function(refine,seed) {
  n = 6 * refine
  m = 5 * refine
  
  N = 5*n
  M = m
  
  p = as.matrix(expand.grid(x = 1:N-1, y = 1:M-1))
  A = matrix(c(1/n,0,1/(2*m),1/m),2,2)
  pA = p %*% A
  circ = function(n) { x = 1:n-1; ifelse(x > n/2, x-n, x)/n }
  f = as.matrix(expand.grid(x = circ(N), y = circ(M)))
  fA = f %*% solve(t(A))
  
  sel = rowSums(abs(round(fA) - fA) > 1e-5) == 0
  
  k = sqrt(rowSums(fA^2))
  
  K = matrix(0,N,M)
  K[sel] = (scale/(k^alpha))[sel]
  K[1 ] = 0
  
  corr = corr.profile(1/k)
  a = atan(corr)
  
  RN = ordered_rnorm_mat(N,M,3,seed=seed)
  f1 = fft(K*(cos(a)*RN[[1]] + sin(a)*RN[[2]]),inverse=TRUE)
  f2 = fft(K*(sin(a)*RN[[1]] + cos(a)*RN[[2]]),inverse=TRUE)
  diffvar = sum(2 * K^2 * (cos(a)-sin(a))^2)
  if (is.na(gap)) gap = -qnorm(closed,mean=0,sd=sqrt(diffvar))
  f1 = Re(f1) + gap/2
  f2 = Re(f2) - gap/2
    
  list(f1 = f1,
       f2 = f2)
}
    

seam.geom = function(refine=1,
                     generator,
                     seed=0) {
  n = 6 * refine
  m = 5 * refine
  
  N = 5*n
  M = m
  
  R = matrix(0,N,M)
  X = row(R)-1
  Y = col(R)-1
  B = Y/m + X/m / 2
  A = X/n
  
  ret = generator(refine,seed)
  f1 = ret$f1
  f2 = ret$f2
  ret$f1 = NULL
  ret$f2 = NULL
  
  A = A[1:(2*n),1:m]
  B = B[1:(2*n),1:m]
  f1 = f1[1:(2*n),1:m]
  f2 = f2[1:(2*n),1:m]
  A = cbind(A,A,A,A)
  B = cbind(B-2,B-1,B,B+1)
  f1 = cbind(f1,f1,f1,f1)
  f2 = cbind(f2,f2,f2,f2)
  
  #l = quantile(f1-f2,0.1)
  #f1 = f1 - l/2
  #f2 = f2 + l/2
  
  P = data.frame(x = as.vector(A),y = as.vector(B))
  P$f1 = as.vector(f1)
  P$f2 = as.vector(f2)
  I = 1:length(A)
  dim(I) = dim(A)
  i = rbind(cbind(
    as.vector(I[-nrow(I),-ncol(I)]),
    as.vector(I[-1,-ncol(I)]),
    as.vector(I[-nrow(I),-1])
  ),cbind(
    as.vector(I[-1,-1]),
    as.vector(I[-nrow(I),-1]),
    as.vector(I[-1,-ncol(I)])
  ))
  
  sel = P$y < 2 & P$y >= 0
  
  ni = rep(0,nrow(P))
  ni[sel] = 1:sum(sel)
  i[] = ni[i]
  i = i[rowSums(i == 0) == 0,]
  P = P[sel,]
  
  #iv = as.vector(t(i))
  #clear3d()
  #triangles3d(P1[iv,],col=2)
  #triangles3d(P2[iv,],col=3)
  
  sel = P$f1 < P$f2
  P$fm = (P$f1 + P$f2)/2
  P$f1[sel] = P$fm[sel]
  P$f2[sel] = P$fm[sel]
  
  sel = sel[i]
  dim(sel) = dim(i)
  i = i[rowSums(sel) != 3,]
  sel = rep(FALSE,nrow(P))
  sel[i[,1]] = TRUE
  sel[i[,2]] = TRUE
  sel[i[,3]] = TRUE
  
  ni = rep(0,nrow(P))
  ni[sel] = 1:sum(sel)
  i[] = ni[i]
  i = i[rowSums(i == 0) == 0,]
  P = P[sel,]
  P$x = P$x - 0.5
  P$y = P$y - 0.5
  P$h = P$f1 - P$f2
  
  bonds2 = range(P$f1,P$f2)
  cat("Final Bonds:",bonds2[1],bonds2[2],"\n")
  
  ret$points=P
  ret$triangles=i
  return(ret)
}


seam3d = function(obj,type=c("top","bottom"),top="top" %in% type,bottom="bottom" %in% type,middle="middle" %in% type) {
  iv = as.vector(t(obj$triangles))
  clear3d()
  plot3d(diag(3))
  clear3d()
  if (top) triangles3d(obj$points$f1[iv],obj$points$x[iv],obj$points$y[iv],col=2)
  if (bottom) triangles3d(obj$points$f2[iv],obj$points$x[iv],obj$points$y[iv],col=3)
  if (middle) triangles3d(obj$points$fm[iv],obj$points$x[iv],obj$points$y[iv],col=4)
}

save.msh = function(obj, filename,type=c("top","bottom"),top="top" %in% type,bottom="bottom" %in% type,middle="middle" %in% type) {
  points = NULL
  triangles = NULL
  i = 0
  n = nrow(obj$points)
  if (top) {
    points = rbind(points, obj$points[,c("f1","x","y")])
    triangles = rbind(triangles, obj$triangles + i*n)
    i = i + 1
  }
  if (bottom) {
    points = rbind(points, obj$points[,c("f2","x","y")])
    triangles = rbind(triangles, obj$triangles + i*n)
    i = i + 1
  }
  if (middle) {
    points = rbind(points, obj$points[,c("fm","x","y")])
    triangles = rbind(triangles, obj$triangles + i*n)
    i = i + 1
  }
  f = file(filename,"w")
  cat("Triangles\n",file=f)
  cat("3D-Nodes ",nrow(obj$points),"\n",sep="",file=f)
  i = 1:nrow(points)-1
  write.table(cbind(i,i,0,points),file=f,row.names=FALSE,col.names=FALSE,sep="\t")
  cat("\nTri3 ",nrow(triangles),"\n",file=f,sep="")
  i = 1:nrow(triangles)-1
  write.table(cbind(i,0,triangles-1),file=f,row.names=FALSE,col.names=FALSE,sep="\t")
  close(f)
}

seam.cut = function(obj, eps = 1e-9){
  snap = function(x) ifelse(x > -eps & x < eps, 0, ifelse(x > 1-eps & x < 1+eps, 1, x))
  i = obj$triangles
  sel = obj$points$x >= -eps & obj$points$x <= 1+eps & obj$points$y >= -eps & obj$points$y <= 1+eps
  sel = sel[i]
  dim(sel) = dim(i)
  tocut = (! sel[,1]) & sel[,2] & sel[,3]
  obj$points[i[tocut,1],] = (obj$points[i[tocut,1],] + obj$points[i[tocut,3],])/2
  
  sel = obj$points$x >= -eps & obj$points$x <= 1+eps & obj$points$y >= -eps & obj$points$y <= 1+eps
  ni = rep(0,nrow(obj$points))
  ni[sel] = 1:sum(sel)
  i[] = ni[i]
  i = i[rowSums(i == 0) == 0,]
  ret = obj
  ret$points = obj$points[sel,]
  ret$points$x = snap(ret$points$x)
  ret$points$y = snap(ret$points$y)
  ret$triangles = i
  ret
}

seam.volume = function(obj) {
  v = obj$points[,c("x","y")]
  v1 = v[obj$triangles[,2],] - v[obj$triangles[,1],]
  v2 = v[obj$triangles[,3],] - v[obj$triangles[,1],]
  a = abs(1/2*(v1[,1]*v2[,2] - v1[,2]*v2[,1]))
  h = 1/3*(obj$points$h[obj$triangles[,1]] + obj$points$h[obj$triangles[,2]] + obj$points$h[obj$triangles[,3]])
  sum(a*h)
}
