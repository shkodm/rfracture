
phip = function(x,a,b){
  ifelse(x<a,0,ifelse(x>b,1,(log(x)-log(a))/(log(b)-log(a))))
}

dhip = function(x,a,b){
  ifelse(x<a,0,ifelse(x>b,0,1/x*1/(log(b)-log(a))))
}

qhip = function(q,a,b){
  ifelse(q<0,a,ifelse(q>1,b,(b/a)^q*a))
}

rhip = function(n,a,b){
  qhip(runif(n),a,b)
}


