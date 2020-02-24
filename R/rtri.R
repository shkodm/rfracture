rtri = function(n, h1,h2,h3) {
  a = runif(n)
  b = runif(n)
  sel = a+b > 1
  a = ifelse(sel,1-a,a)
  b = ifelse(sel,1-b,b)
  c = 1-a-b
  h = runif(n,min = 0,max=h1+h2+h3)
  ht = a * h1 + b*h2 + c*h3
  sel = h > ht
  h[sel] = h[sel] - ht[sel]
  pom = a[sel]
  a[sel] = b[sel]
  b[sel] = c[sel]
  c[sel] = pom
  ht = a * h1 + b*h2 + c*h3
  sel = h > ht
  h[sel] = h[sel] - ht[sel]
  pom = a[sel]
  a[sel] = b[sel]
  b[sel] = c[sel]
  c[sel] = pom
  ht = a * h1 + b*h2 + c*h3
  sel = h > ht
  if (any(sel)) stop("Something it wrong")
  data.frame(w1=a,w2=b,w3=c,ht=ht,h=h,w4=h/ht)
}

#ret = rtri(20000,-1,2,2)
#plot(ret$w1,ret$w2,asp=1,xlim=c(0,1),ylim=c(0,1))
#library(rgl)
#plot3d(ret$w1,ret$w2,ret$h,xlim=c(0,1),ylim=c(0,1))
