
#########################################################################
#                                                                       #
#       Function1:                                                      #
#       Create matrix by repeating rows                                 #
#                                                                       #
#########################################################################
rep.row <- function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#########################################################################
#                                                                       #
#       Function2:                                                      #
#       Bootstrap control limit for phase0 T2                           #
#                                                                       #
#########################################################################
BootstrapUCL0 <- function(x0, alpha, iterations){
  N0 <- sapply(x0, nrow) %>% sum # number of observations in phase0
  UCL <- vector("numeric")
  for (i in 1:iterations) {
    x0PooledSampled <- x0 %>% do.call(what=rbind) %>% # pool observations in phase0
      .[sample(N0, N0, replace=T),] # Sample x0 matrix by bootstrapping it's rows
    
    UCL[i] <- x0PooledSampled %>% apply(FUN=mean, MARGIN=2) %>%
      rep.row(N0) %>% -(x0PooledSampled) %>%
      apply(MARGIN=1, function(x) {x %*% ginv(cov(x0PooledSampled)) %*% x}) %>% # T2 in phase0
      quantile(probs=1-alpha)
  }
  return(mean(UCL)) # Mean of 95th percentile of UCL's
}

#########################################################################
#                                                                       #
#       Function3:                                                      #
#       Initialize the decomposition matrix                             #
#                                                                       #
#########################################################################
# Taken from ellipseChart function, qcc package
init.decomp <- function(p){
  k=0
  for(i in 1:p) {
    k <- k+factorial(p)/(factorial(i)*factorial(p-i)) # total number of permutations 
  }
  
  for (i in 1:k) {
    v=1
    q <- matrix(0, k, p + 2)  # matrix of all the permutations
    for(i in 1:p) { 
      a <- t(combn(p, i))
      for(l in 1:nrow(a)) {
        for(j in 1:ncol(a)) {
          q[v, j+2] <- a[l, j]
        }
        v = v + 1
      }
    }
  }
  return(q)
}









