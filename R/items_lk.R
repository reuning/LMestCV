items_lk <- function(S,R,yv,piv,Pi,Psi, Ug){

# Preliminaries
  	sS = dim(S)
  	ns = sS[1]
  	TT = sS[2]
  	k = max(Ug)
  	if(length(sS)==2) r = 1 else r = sS[3]
  	if(r==1){
  		if(is.matrix(S)) S = array(S,c(dim(S),1))
  		if(is.matrix(R)) R = array(R,c(dim(R),1))
  	}
  	miss = !is.null(R)
# Compute log-likelihood
  	Phi = array(NA,c(ns,TT))
	L = array(0,c(ns,TT))
	if(miss) {
	  for(j in 1:r){
	    tmp = (Psi[S[,1,j]+1,,j]*R[,1,j]+(1-R[,1,j]))
	    Phi[,1] = tmp[cbind(1:ns, Ug[,1])]
	  } 
	}	else {
	  for(j in 1:r) Phi[,,1] = Phi[,1]*Psi[S[,1,j]+1,Ug[,j],j]
	}
  	L[,1] = Phi[,1]
  	for(t in 2:TT){
  		if(miss) {
  		  for(j in 1:r){
  		    tmp = (Psi[S[,t,j]+1,,j]*R[,t,j]+(1-R[,t,j]))
  		    Phi[,t]= tmp[cbind(1:ns, Ug[,t])]
  		  } 
  		} else{
  		  for(j in 1:r) Phi[,t] = Psi[S[,t,j]+1,,j]
  		} 
   		L[,t] = Phi[,t]
  	}
  	if(ns==1) {
  	  pv = sum(L)
  	} else{
  	  pv = apply(L, 1, prod)
  	} 
  	lk = sum(yv*log(pv))
  	out = list(lk=lk,Phi=Phi,L=L,pv=pv)
}
