
#' Generate a matrix of random complex numbers in a consistent order
#' 
#' @param N number of rows
#' @param M number of columns
#' @param k number of independent matrices
#' @param seed random seed
#' @param length_one if TRUE, return complex numbers of length 1
#' @example 
#' ordered_rnorm_mat(5,5,1,seed=123)
#' @export
ordered_rnorm_mat = function(N, M, k=1, seed, length_one=FALSE) {
  if (! missing(seed)) set.seed(seed)
  RN = t(matrix(rnorm(N*M*k*2),k*2,N*M))
  R = matrix(0+0i,N,M)
  circ = function(x,n) ifelse(x > n/2,x-n,x)
  fpZ = expand.grid(i=circ(1:N-1,N)/N,j=circ(1:M-1,M)/M)
  o = order(pmax(abs(fpZ$i),abs(fpZ$j)),atan2(fpZ$i,fpZ$j))
  lapply(1:k-1, function(i) {
    z = RN[,2*i+1] + 1i*RN[,2*i+2]
    R[o] = z
    if (length_one) {
      I = matrix(1:length(R),nrow(R),ncol(R))
      I2 = I[c(1,nrow(R):2),c(1,ncol(R):2)]
      R2 = Conj(R[c(1,nrow(R):2),c(1,ncol(R):2)])
      R[] = ifelse(I == I2, 0, (R + R2)/Mod(R+R2))
    }
    R
  })
}

fracture.geom.my = function(totalScale = 0.2,
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

#' Correlation profile with linear cut-off
#' 
#' @param ML Length of correlation cut-off
#' @param TL Width of linear cut-off fade
#' @param MinMF Minimal correlation
#' @param MaxMF Maximal correlation
#' @references
#' Steven R. Ogilvie, Evgeny Isakov, Paul W.J. Glover,
#' Fluid flow through rough fractures in rocks. II: A new matching model for rough rock fractures,
#' Earth and Planetary Science Letters, Volume 241, Issues 3–4, 2006, https://doi.org/10.1016/j.epsl.2005.11.041.
#' 
#' @export
ogilvie.corr.profile = function(ML=0.5, TL=0, MinMF=0, MaxMF=1) function(lambda) {
  ifelse(lambda >= ML+TL/2, 1, ifelse(lambda <= ML-TL/2, 0, (ML-TL/2-lambda)*(lambda-(ML+3*TL/2))/(TL*TL) ))*(MaxMF-MinMF)+MinMF
}

# spectral density based on Brown: Simple mathematical model of a rough fracture
fracture.geom.brown = function(scale=1, alpha=2, corr.profile, closed=0.1, gap=NA, beta.make = FALSE, beta.shape = 4)  function(refine,seed) {
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
  f1 = Re(f1)
  f2 = Re(f2)
  diff_sd = sqrt(sum(2 * K^2 * (cos(a)-sin(a))^2))
  sum_sd = sqrt(sum(2 * K^2 * (cos(a)+sin(a))^2))
  if (beta.make) {
    fsum = f1 + f2
    fdiff = f1 - f2
    beta_sd = sqrt(1/(4*(2*shape+1)))
    cat("Centerline scale:",1/beta_sd*sum_sd/2,"\n")
    fsum = (qbeta(pnorm(fsum,sd=sum_sd),shape1=shape,shape2=shape)-0.5)/beta_sd*sum_sd
    cat("Gap scale:",1/beta_sd*diff_sd,"\n")
    fdiff = (qbeta(pnorm(fdiff,sd=diff_sd),shape1=shape,shape2=shape)-qbeta(closed,shape1=shape,shape2=shape))/beta_sd*diff_sd
    f1 = (fsum + fdiff)/2
    f2 = (fsum - fdiff)/2
    max_fdiff = (1-qbeta(closed,shape1=shape,shape2=shape))/beta_sd*diff_sd
    cat("Max gap:", max_fdiff,"\n")
    cat("Bonds from beta:",(0.5/beta_sd*sum_sd + max_fdiff)/2,-(0.5/beta_sd*sum_sd + max_fdiff)/2,"\n")
  } else {
    if (is.na(gap)) gap = -qnorm(closed,mean=0,sd=diff_sd)
    f1 = f1 + gap/2
    f2 = f2 - gap/2
  }
  list(f1 = f1,
       f2 = f2)
}


plot.fracture_field = function(obj, field="f1", col.palette=c("black","red","yellow"), pch=16, cex=1, asp=1, ...){
  col = as.vector(obj[[field]])
  col = (col-min(col))/(max(col)-min(col))
  col = colorRamp(col.palette)(col)
  col = rgb(col,max=255)
  plot(obj$points[,1], obj$points[,2], col=col,pch=pch,cex=cex,asp=asp,...)
  abline(h=-5:5); abline(v=-5:5)
  arrows(0,0,obj$mat[,1],obj$mat[,2],angle = 15)
}


exp.spectrum = function(scale=1, alpha=2) function(k) scale/(k^alpha)

fracture.geom.fields = function(
    N, M,
    mat = diag(c(1,1)),
    A = diag(c(1/N,1/M)) %*% mat,
    spectrum = exp.spectrum(),
    corr.profile = function(k) 0,
    closed = 0.1,
    gap = NA,
    seed,
    length_one = FALSE,
    bonds,
    cut=TRUE,
    widen=0, widen_grad=1) {
  p = as.matrix(expand.grid(x = 1:N-1, y = 1:M-1))
  pA = p %*% A
  circ = function(n) { x = 1:n-1; ifelse(x > n/2, x-n, x)/n }
  f = as.matrix(expand.grid(x = circ(N), y = circ(M)))
  fA = f %*% solve(t(A))
  
  sel = rowSums(abs(round(fA) - fA) > 1e-6) == 0
  Scales = data.frame(
    i  = (1:nrow(fA))[sel],
    kx = as.integer(fA[sel,1]),
    ky = as.integer(fA[sel,2]))
  Scales$k = sqrt(Scales$kx^2 + Scales$ky^2)
  Scales = Scales[Scales$k != 0,]
  Scales$Power = spectrum(Scales$k)
  
  MaxK = max(Scales$kx, Scales$ky, -Scales$kx, -Scales$ky)
  
  RN = ordered_rnorm_mat(MaxK*2+1,MaxK*2+1,3,seed=seed,length_one = length_one)
  Scales$RNIndex = ifelse(Scales$kx<0, MaxK*2+1 +Scales$kx, Scales$kx) + ifelse(Scales$ky<0, MaxK*2+1 +Scales$ky, Scales$ky)*(MaxK*2+1) + 1
  Scales$corr = corr.profile(1/Scales$k)
  Scales$corr.angle = atan(Scales$corr)
  Scales$PowerR1F1 = Scales$Power * cos(Scales$corr.angle)
  Scales$PowerR2F1 = Scales$Power * sin(Scales$corr.angle)
  Scales$PowerR1F2 = Scales$Power * sin(Scales$corr.angle)
  Scales$PowerR2F2 = Scales$Power * cos(Scales$corr.angle)
  
  K = matrix(0,N,M)
  K[Scales$i] = Scales$PowerR1F1 * RN[[1]][Scales$RNIndex] + Scales$PowerR2F1 * RN[[2]][Scales$RNIndex]
  f1 = Re(fft(K,inverse=TRUE))
  K[Scales$i] = Scales$PowerR1F2 * RN[[1]][Scales$RNIndex] + Scales$PowerR2F2 * RN[[2]][Scales$RNIndex]
  f2 = Re(fft(K,inverse=TRUE))
  
  Scales$PowerR1Diff = Scales$PowerR1F1 - Scales$PowerR1F2
  Scales$PowerR2Diff = Scales$PowerR2F1 - Scales$PowerR2F2
  diff_sd = sqrt(sum(Scales$PowerR1Diff^2 + Scales$PowerR2Diff^2))
  f1_sd = sqrt(sum(Scales$PowerR1F1^2 + Scales$PowerR2F1^2))
  f2_sd = sqrt(sum(Scales$PowerR1F2^2 + Scales$PowerR2F2^2))
  f1_mean = 0
  f2_mean = 0
  if (is.na(gap)) gap = -qnorm(closed,mean=0,sd=diff_sd)
  f1 = f1 + gap/2
  f1_mean = f1_mean + gap/2
  f2 = f2 - gap/2
  f2_mean = f2_mean - gap/2
  
  sel = f1 < f2
  fm = (f1 + f2)/2
  f1[sel] = fm[sel]
  f2[sel] = fm[sel]
  h2 = f1 - fm
  widen_fac = (pnorm(h2,sd=1/(sqrt(2*pi)*widen_grad/(2*widen)))*2-1)*widen
  f1 = f1 + widen_fac
  f2 = f2 - widen_fac
  if (!missing(bonds)) { if (length(bonds) == 2) {
    shift = mean(bonds)
    f1 = f1 + shift
    f1_mean = f1_mean + shift
    f2 = f2 + shift
    f2_mean = f2_mean + shift
    Pcut = c(pnorm(bonds,mean=f1_mean,sd=f1_sd),pnorm(bonds,mean=f2_mean,sd=f2_sd))
    Pcut = max(Pcut[1],1-Pcut[2],Pcut[3],1-Pcut[4])
    cat("Probability of being out of bonds:", Pcut,"\n")
    if (cut) {
      if (Pcut > 1e-2) warning("Probability of being out of bonds is higher then 1%")
      f1[] = ifelse(f1<bonds[1],bonds[1],ifelse(f1>bonds[2],bonds[2],f1))
      f2[] = ifelse(f2<bonds[1],bonds[1],ifelse(f2>bonds[2],bonds[2],f2))
    }
  } else stop("bonds should be of length 2")}
  ret = list(f1 = f1, f2 = f2, points = pA, scales = Scales, mat=mat, A=A, f1_mean=f1_mean, f1_sd=f1_sd, f2_mean=f2_mean, f2_sd=f2_sd)
  class(ret) = "fracture_field"
  ret
}


fracture.geom = function(refine=1, ...) {
  n = 6 * refine
  m = 5 * refine
  
  N = 5*n
  M = m

  ret = fracture.geom.fields(N, M, mat = matrix(c( 5, 0, 3, 1),2,2), ...)

  f1 = ret$f1
  f2 = ret$f2
  ret$f1 = NULL
  ret$f2 = NULL
  A = ret$points[,1]
  dim(A) = dim(f1)
  B = ret$points[,2]
  dim(B) = dim(f1)
  
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
  

  P$x = P$x - 0.5
  P$y = P$y - 0.5
  P$h = P$f1 - P$f2
  
  bonds2 = range(P$f1,P$f2)
  cat("Final Bonds:",bonds2[1],bonds2[2],"\n")
  
  ret$points=P
  ret$triangles=i
  return(ret)
}


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

#' Makes a 3d plot of the fracture
#'
#' @export
#' @import rgl
fracture3d = function(obj,type=c("top","bottom"),top="top" %in% type,bottom="bottom" %in% type,middle="middle" %in% type,col=c(2,3,4), add=FALSE) {
  if (length(col) == 1) col = rep(col,3)
  iv = as.vector(t(obj$triangles))
  if (!add) {
    clear3d()
    plot3d(diag(3))
    clear3d()
  }
  if (top) triangles3d(obj$points$f1[iv],obj$points$x[iv],obj$points$y[iv],col=col[1])
  if (bottom) triangles3d(obj$points$f2[iv],obj$points$x[iv],obj$points$y[iv],col=col[2])
  if (middle) triangles3d((obj$points$f1+obj$points$f2)[iv]/2,obj$points$x[iv],obj$points$y[iv],col=col[3])
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

fracture.cut = function(obj, eps = 1e-9){
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

fracture.volume = function(obj) {
  v = obj$points[,c("x","y")]
  v1 = v[obj$triangles[,2],] - v[obj$triangles[,1],]
  v2 = v[obj$triangles[,3],] - v[obj$triangles[,1],]
  a = abs(1/2*(v1[,1]*v2[,2] - v1[,2]*v2[,1]))
  h = 1/3*(obj$points$h[obj$triangles[,1]] + obj$points$h[obj$triangles[,2]] + obj$points$h[obj$triangles[,3]])
  sum(a*h)
}


border3d = function(obj, f1, f2, add=FALSE, ...) {
  edges = rbind(obj$triangles[,1:2],obj$triangles[,2:3],obj$triangles[,c(1,3)])
  sel = edges[,1] > edges[,2]
  edges[sel,] = edges[sel,c(2,1)]
  head(edges)
  edges = edges[order(edges[,1],edges[,2]),]
  a = !duplicated(edges)
  head(edges[a,])
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
    plot3d(diag(3))
    clear3d()
  }
  triangles3d(rbind(t1,t2),...)
}
