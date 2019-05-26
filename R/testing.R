


r <- 6 #### number of items
b <- 7 #### number of categorties
k <- 6 #### number of classes
n <- 1000 ### number of responses
TT <- 3 ### number of transitions


library(gtools)

piv <- rep(1/k, k)
Pi <- array((1/5)*(1/k), dim=c(k, k, TT))
Pi[,,1] <- 0
for(ii in 2:TT) diag(Pi[,,ii]) <- ((4*k+1)/5)*(1/k)

Psi <- c(replicate(k*r, rdirichlet(1, alpha=sample(c(rep(.5, b-1), 10), size=b))))
Psi <- array(Psi, c(b, k, r))
# Psi <- apply(Psi, c(2,3), function(x) exp(x)/sum(exp(x)))


sample_Y <- draw_lm_basic_eff(piv, Pi, Psi, n)

sample_Y = as.data.table(sample_Y)
sample_Y = sample_Y[,.N, by=names(sample_Y)]
S_boot = as.matrix(sample_Y)[,-ncol(sample_Y)]
yv_boot = sample_Y$N
S_boot = array(t(S_boot), c(r, TT, length(yv_boot)))
S_boot = aperm(S_boot)
if (r == 1)  S_boot = S_boot[, , 1]


check <- search.model.LM(version="basic", kv=1:10, S=S_boot, yv=yv_boot, mod=0, out_se=F)

plot(check$aicv)
which.min(check$aicv)

ll_sim <- matrix(NA, 10, r)
aic_sim <- matrix(NA, 10, r)
bic_sim <- matrix(NA, 10, r)

for(kk in 2:10){
  
  tmp_sim <- item_heldout(S_boot, yv_boot,k=kk, mod=0)
  ll_sim[kk,] <- tmp_sim$lk
  aic_sim[kk,] <- tmp_sim$aic
  bic_sim[kk,] <- tmp_sim$bic
  
}

plot(rowSums(aic_sim))
which.min(rowSums(aic_sim))

plot(rowSums(bic_sim))
which.min(rowSums(bic_sim))
