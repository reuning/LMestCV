#' checks heldout
#'@export
item_heldout <- function(S, yv, k, mod=1, out_se=F, full=F){
  
  S_full <- S
  items <- dim(S_full)[3]
  
  ll_out <- matrix(NA, length(items))
  aic_out <- matrix(NA, length(items))
  bic_out <- matrix(NA, length(items))
  
  cat("------------------------------------------------------------\n")
  
  for(jj in 1:items){
    
    
    S <- S_full[,,-jj]
    sink("/dev/null")
    if(full){
      mod_full <- find_lm_basic(S=S, yv=yv, k=k, mod=mod, out_se=F)
      mod_full <- mod_full$out.single
    } else {
      mod_full <- est_lm_basic(S=S, yv=yv, k=k, mod=mod, out_se=F)
      
    }
    
    S <- S_full[,,jj]
    if(full){
      mod_tmp <- est_lm_cont(S=S, yv=yv, k=k, mod=mod, out_se=F, Pi=mod_full$Pi, 
                             piv = mod_full$piv)
      mod_tmp <- mod_tmp$out.single
      
    } else {
      mod_tmp <- est_lm_cont(S=S, yv=yv, k=k, mod=mod, out_se=F, Pi=mod_full$Pi, 
                             piv = mod_full$piv)
    }
      

    sink()
    
    cat("item:", jj, "loglik:",mod_tmp$lk, "AIC:", mod_tmp$aic, "BIC:", mod_tmp$bic, '\n' )
    ll_out[jj] <- mod_tmp$lk
    aic_out[jj] <- mod_tmp$aic
    bic_out[jj] <- mod_tmp$bic
    
  }
  cat("------------------------------------------------------------\n")
  
  out <- list("lk"=ll_out, "aic"=aic_out, "bic"=bic_out)
  return(out)
}