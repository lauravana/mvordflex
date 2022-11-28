## Issue: the entries of the
set.seed(123)
PL.Q <- function(P, L) {
  A <- P %*% L
  QR <- base::qr(A)
  R  <- base::qr.R(QR)
  Q <- A %*% solve(R)
}
for (i in 1:100){
  m <- 3
  ## Following Roy 2019
  ## unrestricted parameters:
  l <- rnorm(m * (m - 1)/2) ## cholesky of V positive definite
  d <- rnorm(m)             ## cholesky of V positive definite
  s <- rnorm(m * (m - 1)/2)        ## Q orthogonal
  ## Transformations
  ### For V - positive definite
  L <- diag(m)
  #L[lower.tri(L)] <- l
  #D <- diag(exp(d))
  #V <- L %*% tcrossprod(D, L)
  # V <- L %*% (t(L) * exp(d))
  #L[lower.tri(L)] <- l
  #diag(L) <- exp(d)
  #V <- tcrossprod(L)
  ## log param
  logV <- matrix(0, nrow = m, ncol = m)
  logV[upper.tri(logV)] <- l#c(l, d)
  logV <- logV + t(logV)
  diag(logV) <- d
  e <- eigen(logV)
  U <- e$vectors
  loglambda <- e$values
  V <- U %*% diag(exp(loglambda)) %*%t(U)
  ### For Q - orthogonal
  #S <- matrix(0, nrow = m, ncol = m)
  S <- diag(m)
  S[lower.tri(S)] <- s
  #S <- S - t(S)
  # S is skeq symmetric: S^T = -S
  #  Q  = (diag(m) - S) (diag(m) + S) ^ {-1}
  Im <- diag(m)
  theta <- pi * (1 - exp(s))/(1 + exp(s))
  #Q <- Im#PL.Q(Im[, sample(1:m)], S)
  index<-1
  U_transpose <- diag(m)
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      cos1=cos(theta[index]);
      sin1=sin(theta[index]);
      Ui_tmp=cos1*U_transpose[,i]-sin1*U_transpose[,j];
      Uj_tmp=sin1*U_transpose[,i]+cos1*U_transpose[,j];
      U_transpose[,i]=Ui_tmp;
      U_transpose[,j]=Uj_tmp;
      index <- index+1;
    }
    #index<-1;
  }
  Q <- (U_transpose)
  theta <- pi * (1 - exp(s))/(1 + exp(s))
  ind <- combn(m, 2)
  G <- diag(m)
  for (k in ncol(ind)){
    i <- ind[1, k]
    j <- ind[2, k]
    Gk <- diag(m)
    Gk[i, i] <- Gk[j, j] <- cos(theta[k])
    Gk[i, j] <- - sin(theta[k])
    Gk[j, i] <- + sin(theta[k])
    G <- G %*% Gk
  }
  Q<-G
  #Q <- solve(Im - S, Im + S)
  #Q <- (diag(m) - S) %*% solve(diag(m) + S)
  ## make matrix A
  M <- diag(m)
  Vhalf <- U %*% diag(exp(loglambda/2)) %*%t(U)
  Psi <- Vhalf %*% Q %*% amen::mhalf(chol2inv(chol(V + M)))
  #Psi <- chol(V) %*% Q %*% chol(chol2inv(chol(V + M)))

  #Psi <- chol(V) %*% Q %*% chol(chol2inv(chol(V + M)))
  print(eigen(Psi)$val)
  print(all(abs(eigen(Psi)$val) < 1))
}

library("amen")

mhalf(V)
mhalf(Q)

A <- mhalf(V) %*% Q %*% mhalf(chol2inv(chol(V + M)))
A
