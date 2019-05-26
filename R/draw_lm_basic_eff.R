#' Draws samples with a more efficient version of data agregation
#'@export
draw_lm_basic_eff <-
function(piv,Pi,Psi,n){

#        [lk,piv,Pi,Psi,np,aic,bic,lkv] = draw_lm_basic(piv,Pi,Psi,n)
#
# Draw a sample of size n from a Basic Latent Markov model with parameter piv, Pi and Psi

# Preliminaries
  k = length(piv)
  dd = dim(Psi)
  c = apply(Psi, c(2,3), function(x) sum(!is.na(x)))[1,]
  TT = dim(Pi)[3]
  if (length(dd) > 2) 
    r = dd[3]
  else r = 1
  Y = matrix(0, n, TT * r)
  cat("------------|\n")
  cat(" sample unit|\n")
  cat("------------|\n")
  for (i in 1:n) {
    if (i/1000 == floor(i/1000)) 
      cat(sprintf("%11g", i), "\n", sep = " | ")
    u = k + 1 - sum(runif(1) < cumsum(piv))
    ind = 0
    for (j in 1:r) {
      ind = ind + 1
      Y[i, ind] = c[j] - sum(runif(1) < cumsum(Psi[, u, j]), na.rm=T)
    }
    for (t in 2:TT) {
      u = k + 1 - sum(runif(1) < cumsum(Pi[u, , t]))
      for (j in 1:r) {
        ind = ind + 1
        Y[i, ind] = c[j] - sum(runif(1) < cumsum(Psi[, 
                                                     u, j]), na.rm=T)
      }
    }
  }
	cat("------------|\n")
	return(Y)

}