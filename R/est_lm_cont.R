#' contrainsted estimations
#'@export
est_lm_cont <-
  function(S,
           yv,
           k,
           start = 0,
           mod = 0,
           tol = 10 ^ -8,
           maxit = 1000,
           out_se = FALSE,
           piv = NULL,
           Pi = NULL,
           Psi = NULL) {
    # Preliminaries
    check_der = FALSE  # to check derivatives
    n = sum(yv)
    sS = dim(S)
    ns = sS[1]
    TT = sS[2]
    if (min(S, na.rm = T) > 0) {
      cat("|------------------- WARNING -------------------|\n")
      cat("|The first response category must be coded as 0 |\n")
      cat("|-----------------------------------------------|\n")
    }
    
    if (is.data.frame(S)) {
      warning("Data frame not allowed for S")
    }
    
    if (ns != length(yv))
      stop("dimensions mismmissatch between S and yv")
    
    if (length(sS) == 2) {
      r = 1
      if (is.matrix(S))
        S = array(S, c(dim(S), 1))
    } else
      r = sS[3]
    
    miss = any(is.na(S))
    if (miss) {
      cat("Missing data in the dataset, treated as missing at random\n")
      R = 1 * (!is.na(S))
      S[is.na(S)] = 0
    } else{
      R = NULL
    }
    Sv = matrix(S, ns * TT, r)
    if (miss)
      Rv = matrix(R, ns * TT, r)
    bv = apply(Sv, 2, max) ### number of catgories in each 
    b = max(bv) ## maximum number of categories
    m = vector("list", r) 
    for (j in 1:r) {
      m[[j]]$Co = cbind(-diag(bv[j]), diag(bv[j]))
      Maj = cbind(lower.tri(matrix(1, bv[j], bv[j]), diag = TRUE), rep(0, bv[j]))
      m[[j]]$Ma = rbind(Maj, 1 - Maj)
    }
    th = NULL
    sc = NULL
    J = NULL
    
    if (k == 1) {
      piv = 1
      Pi = 1
      P = matrix(NA, b + 1, r) ### Number of Responses in each category 
      for (j in 1:r)
        P[1:(bv[j] + 1), j] = 0
      for (t in 1:TT) {
        for (j in 1:r) {
          for (y in 0:b) {
            ind = which(S[, t, j] == y)
            P[y + 1, j] = P[y + 1, j] + sum(yv[ind])
          }
        }
      }
      if (miss){
        P[1,] = P[1,] - colSums(Rv==0) ####### Drops out all the NAs that were added as 0s  
      }
      Psi = apply(P, 2, function(x) x/sum(x, na.rm=T))
      dimnames(Psi) = list(category = 0:b, item = 1:r)
      pm = rep(1, ns)
      for (t in 1:TT)
        for (j in 1:r)
          pm = pm * Psi[S[, t, j] + 1, j]
      lk = sum(yv * log(pm))
      np = r * b
      aic = -2 * lk + np * 2
      bic = -2 * lk + np * log(n)
      out = list(
        lk = lk,
        piv = piv,
        Pi = Pi,
        Psi = Psi,
        np = np,
        aic = aic,
        bic = bic,
        lkv = NULL,
        J = NULL,
        V = NULL,
        th = NULL,
        sc = NULL,
        call = match.call()
      )
      class(out) = "LMbasic"
      return(out)
    }
    # Starting values
    if (start == 0) {
      P = matrix(NA, b + 1, r) ### Number of Responses in each category 
      E = matrix(NA, b, r)
      for (j in 1:r)
        P[1:(bv[j] + 1), j] = 0
      for (t in 1:TT)
        for (j in 1:r)
          for (y in 0:b) {
            ind = which(S[, t, j] == y)
            P[y + 1, j] = P[y + 1, j] + sum(yv[ind])
          }
      if (miss){
        P[1,] = P[1,] - colSums(Rv==0) ####### Drops out all the NAs that were added as 0s  
      }
      P[P==0] <- min(P[P!=0]) 
      for (j in 1:r){
        E[1:bv[j], j] = m[[j]]$Co %*% log(m[[j]]$Ma %*% P[1:(bv[j] + 1), j])
      }

      Psi = array(NA, c(b + 1, k, r))
      Eta = array(NA, c(b, k, r))
      grid = seq(-k, k, 2 * k / (k - 1))
      for (c in 1:k)
        for (j in 1:r) {
          etac = E[1:bv[j], j] + grid[c]
          Eta[1:bv[j], c, j] = etac
          Psi[1:(bv[j] + 1), c, j] = invglob(etac) ##### global logits i assume to provide some stability to it? 
        }
      # piv = rep(1, k) / k
      # Pi = matrix(1, k, k) + 9 * diag(k)
      # Pi = diag(1 / rowSums(Pi)) %*% Pi
      # 
      # Pi = array(Pi, c(k, k, TT))
      # Pi[, , 1] = 0
    }
    if (start == 1) {
      Psi = array(NA, c(b + 1, k, r))
      for (j in 1:r) {
        Psi[1:(bv[j] + 1), , j] = matrix(runif((bv[j] + 1) * k), bv[j] + 1, k)
        for (c in 1:k)
          Psi[1:(bv[j] + 1), c, j] = Psi[1:(bv[j] + 1), c, j] / sum(Psi[1:(bv[j] +
                                                                             1), c, j])
      }
      # Pi = array(runif(k ^ 2 * TT), c(k, k, TT))
      # for (t in 2:TT)
      #   Pi[, , t] = diag(1 / rowSums(Pi[, , t])) %*% Pi[, , t]
      # Pi[, , 1] = 0
      # piv = runif(k)
      # piv = piv / sum(piv)
    }
    if (start == 2) {
      if (is.null(piv))
        stop("initial value of the initial probabilities (piv) must be given in input")
      if (is.null(Pi))
        stop("initial value of the transition probabilities (Pi) must be given in input")
      if (is.null(Psi))
        stop("initial value of the conditional response probabilities (Psi) must be given in input")
      # piv = piv
      # Pi = Pi
      Psi = Psi
    }
    # Compute log-likelihood
    out = complk(S, R, yv, piv, Pi, Psi, k)
    lk = out$lk
    Phi = out$Phi
    L = out$L
    pv = out$pv
    lk0 = sum(yv * log(yv / n))
    dev = 2 * (lk0 - lk)
    cat(
      "------------|-------------|-------------|-------------|-------------|-------------|-------------|\n"
    )
    
    cat("     mod    |      k      |    start    |     step    |     lk      |    lk-lko   | discrepancy |\n")
    
    cat(
      "------------|-------------|-------------|-------------|-------------|-------------|-------------|\n"
    )
    
    cat(sprintf("%11g", c(mod, k, start, 0, lk)), "\n", sep = " | ")
    it = 0
    lko = lk - 10 ^ 10
    lkv = NULL
    par = c(piv, as.vector(Pi), as.vector(Psi))
    if (any(is.na(par)))
      par = par[-which(is.na(par))]
    paro = par
    # Iterate until convergence
    while ((lk - lko) / abs(lk) > tol & it < maxit) {
      Psi0 = Psi
      piv0 = piv
      Pi0 = Pi
      it = it + 1
      
      # ---- E-step ----
      # Compute V and U
      #time = proc.time()
      V = array(0, c(ns, k, TT))
      U = array(0, c(k, k, TT))
      Yvp = matrix(yv / pv, ns, k)
      M = matrix(1, ns, k)
      V[, , TT] = Yvp * L[, , TT]
      U[, , TT] = (t(L[, , TT - 1]) %*% (Yvp * Phi[, , TT])) * Pi[, , TT]
      if (TT > 2) {
        for (t in seq(TT - 1, 2, -1)) {
          M = (Phi[, , t + 1] * M) %*% t(Pi[, , t + 1])
          
          V[, , t] = Yvp * L[, , t] * M
          U[, , t] = (t(L[, , t - 1]) %*% (Yvp * Phi[, , t] * M)) * Pi[, , t]
        }
      }
      M = (Phi[, , 2] * M) %*% t(Pi[, , 2])
      V[, , 1] = Yvp * L[, , 1] * M
      #print(proc.time()-time)
      
      # If required store parameters
      # ---- M-step ----
      # Update Psi
      Y1 = array(NA, c(b + 1, k, r))
      for (j in 1:r)
        Y1[1:(bv[j] + 1)] = 0
      Vv = matrix(aperm(V, c(1, 3, 2)), ns * TT, k)
      for (j in 1:r)
        for (jb in 0:bv[j]) {
          ind = which(Sv[, j] == jb)
          if (length(ind) == 1) {
            if (miss)
              Y1[jb + 1, , j] = Vv[ind, ] * Rv[ind, j]
            else
              Y1[jb + 1, , j] = Vv[ind, ]
          }
          if (length(ind) > 1) {
            if (miss)
              Y1[jb + 1, , j] = colSums(Vv[ind, ] * Rv[ind, j])
            else
              Y1[jb + 1, , j] = colSums(Vv[ind, ])
          }
          
        }
      for (j in 1:r)
        for (c in 1:k) {
          tmp = Y1[1:(bv[j] + 1), c, j]
          if (any(is.na(tmp)))
            tmp[is.na(tmp)] = 0
          tmp = pmax(tmp / sum(tmp), 10 ^ -10)
          Psi[1:(bv[j] + 1), c, j] = tmp / sum(tmp)
        }
      
      #print(proc.time()-time)
      # Update piv and Pi
      # piv = colSums(V[, , 1]) / n
      # U = pmax(U, 10 ^ -300)
      # if (mod == 0)
      #   for (t in 2:TT)
      #     Pi[, , t] = diag(1 / rowSums(U[, , t])) %*% U[, , t]
      # if (mod == 1) {
      #   Ut = apply(U[, , 2:TT], c(1, 2), sum)
      #   Pi[, , 2:TT] = array(diag(1 / rowSums(Ut)) %*% Ut, c(k, k, TT -
      #                                                          1))
      # }
      # if (mod > 1) {
      #   Ut1 = U[, , 2:mod]
      #   if (length(dim(Ut1)) > 2)
      #     Ut1 = apply(Ut1, c(1, 2), sum)
      #   Ut2 = U[, , (mod + 1):TT]
      #   if (length(dim(Ut2)) > 2)
      #     Ut2 = apply(Ut2, c(1, 2), sum)
      #   Pi[, , 2:mod] = array(diag(1 / rowSums(Ut1, 2)) %*% Ut1, c(k, k, mod -
      #                                                                1))
      #   Pi[, , (mod + 1):TT] = array(diag(1 / rowSums(Ut2, 2)) %*% Ut2, c(k, k, TT -
      #                                                                       mod))
      # }
      #print(proc.time()-time)
      # Compute log-likelihood
      paro = par
      par = c(piv, as.vector(Pi), as.vector(Psi))
      if (any(is.na(par)))
        par = par[-which(is.na(par))]
      lko = lk
      
      out = complk(S, R, yv, piv, Pi, Psi, k)
      lk = out$lk
      Phi = out$Phi
      L = out$L
      pv = out$pv
      if (it / 10 == round(it / 10))
        cat(sprintf("%11g", c(mod, k, start, it, lk, lk - lko, max(abs(
          par - paro
        )))), "\n", sep = " | ")
      lkv = c(lkv, lk)
      #print(proc.time()-time)
    }
    # Compute number of parameters
    np =  k * sum(bv)
    # if (mod == 0)
    #   np = np + (TT - 1) * k * (k - 1)
    # if (mod == 1)
    #   np = np + k * (k - 1)
    # if (mod > 1)
    #   np = np + 2 * k * (k - 1)
    aic = -2 * lk + np * 2
    bic = -2 * lk + np * log(n)
    cat(sprintf("%11g", c(mod, k, start, it, lk, lk - lko, max(abs(
      par - paro
    )))), "\n", sep = " | ")
    # adjust output
    if (any(yv != 1))
      V = V / yv
    
    lk = as.vector(lk)
    dimnames(Pi) = list(state = 1:k,
                        state = 1:k,
                        time = 1:TT)
    dimnames(Psi) = list(category = 0:b,
                         state = 1:k,
                         item = 1:r)
    
    out = list(
      lk = lk,
      piv = piv,
      Pi = Pi,
      Psi = Psi,
      np = np,
      aic = aic,
      bic = bic,
      lkv = lkv,
      V = V,
      call = match.call()
    )
    
    cat(
      "------------|-------------|-------------|-------------|-------------|-------------|-------------|\n"
    )
    
    class(out) = "LMbasic"
    return(out)
  }
