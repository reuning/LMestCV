find_lm_basic <- function(S,yv,k,mod=0,nrep=2, tol1 = 10^-5, tol2 = 10^-10, 
                          maxit=1000,out_se=FALSE,piv=NULL,Pi=NULL,Psi=NULL){
  

  cat("***************************************************************************\n")
  out = est_lm_basic(S,yv,k=k,start=0,tol=tol1, mod=0)

  lktrace = out$lk
  lkv = out$lk
  aicv = out$aic
  bicv = out$bic
  cat("lktrace = ",sort(lktrace),"\n")
  cat("lk = ",lkv,"\n")
  cat("aic = ",aicv,"\n")
  cat("bic = ",bicv,"\n")
  
  
  if(k>1){ 
    if(nrep==0){
      cat("***************************************************************************\n")
      cat(c(k,1),"\n")
      outh = est_lm_basic(S,yv,k=k,start=1,tol=tol1, mod=mod)
      
      lktrace = c(lktrace,outh$lk)
      if(outh$lk>out$lk) out = outh	
      
      lkv = out$lk
      aicv = out$aic
      bicv = out$bic
      
      cat("lktrace = ",sort(lktrace),"\n")
      cat("lk = ",lkv,"\n")
      cat("aic = ",aicv,"\n")
      cat("bic = ",bicv,"\n")
    }else{
      for(h in 1:(nrep*(k-1))){
        cat("***************************************************************************\n")
        cat(c(k,h),"\n")
        outh = est_lm_basic(S,yv,k=k,start=1,tol=tol1, mod=mod)

        lktrace = c(lktrace,outh$lk)
        if(outh$lk>out$lk) out = outh	

        lkv = out$lk
        aicv = out$aic
        bicv = out$bic
        
        cat("lktrace = ",sort(lktrace),"\n")
        cat("lk = ",lkv,"\n")
        cat("aic = ",aicv,"\n")
        cat("bic = ",bicv,"\n")
      }
    }
    outn = est_lm_basic(S,yv, mod=mod,k=k,start=2,tol=tol2,piv=out$piv,Pi=out$Pi,Psi=out$Psi,out_se=out_se)
    

    lktrace = c(lktrace,outn$lk)
    out = outn		
    out$lktrace = lktrace
    lkv = out$lk
    aicv = out$aic
    bicv = out$bic
    
  }	
  
  out = list(out.single=out,aicv=aicv,bicv=bicv,lkv=lkv,kv=k,call=match.call())
  return(out)
  
}