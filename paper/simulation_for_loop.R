### Function for creating the (JxTT) x (JxTT) correlation matrix
#' @param blocks a list of matrices
#'               Sigma, Psi * Sigma, Psi^2 * Sigma, ...,Psi^(T-1) * Sigma
toeplitz.block <- function(blocks) {
  l <- length(blocks)
  m.str <- toeplitz(1:l)
  m.str[lower.tri(m.str)] <- 0
  res <- lapply(1:l,function(k) {
    res <- matrix(0, ncol = ncol(m.str), nrow = nrow(m.str))
    res[m.str == k] <- 1
    res %x% blocks[[k]]
  })
  res <- Reduce("+",res)
  res[lower.tri(res)] <-  t(res)[lower.tri(res)]
  res
}
### Simulate one dataset
#' @param J number of ordinal responses
#' @param I number of subjects or firms
#' @param TT number of time points
#' @param P number of covariates to be generated from a standard normal (we assume them constant over the responses)
#' @param beta.matrix a TTxP matrix of coefficients which are allowed to vary over time
#' @param theta a list of length J where each element contains a vector of threshold parameters for response j
#' @param Sigma a JxJ matrix which corresponds to the lag 0 correlation among the J responses
#' @param Psi a JxJ  matrix which in the simplest case is diagonal containing the autoregressive coefficients for each rater.
#'        Otherwise it could be a matrix corresponding to a VAR type of model.
#' @param missings logical indicating whether missing values should be generated at random
#' @param perc.miss percentage of missing values which is constant over all responses

simulate_3dim_1rep_for_loop <- function(J, I, TT, P,
                                        beta.matrix, theta, Sigma, Psi,
                                        missings, perc.miss) {

  firm_id <- rep(1:I, each = TT)
  year_id <- rep(1:TT, I)
  x <- matrix(rnorm(TT * I * P), nrow = TT * I) ## number of covariates including intercept
  X <- array(dim=c(TT * I, P, J))
  for (j in seq_len(J)) X[, , j] <- x
  Xbeta <- sapply(1:J, function(j) {
    xtmat <- model.matrix(~ 0 + factor(year_id):X[,,j])
    xtmat %*% c(beta.matrix)
  })
  # ar_blocks <- vector("list", TT)
  # for (tt in 0:(TT - 1)) {
  #   s_ar <- if (tt == 0) diag(J) else s_ar %*% Psi1
  #   ar_blocks[[tt + 1]] <- s_ar %*% Sigma0
  # }
  #
  # S <- toeplitz.block(ar_blocks)
  ## Simulate errors
  # epsilon <- mvnfast::rmvn(I, rep(0, J * TT), S)
  # errors <- do.call("rbind", lapply(1:nrow(epsilon), function(i)
  #  matrix(epsilon[i, ], byrow = TRUE, ncol = J)))

  epsilon <- lapply(1:I, function(i) {
    tmax <- max(TT, 100)
    u <- MASS::mvrnorm(tmax + 2, mu = rep(0, J), Sigma = Sigma)
    eps <- matrix(nrow = tmax, ncol = J)
    eps0 <-  Psi %*% u[1, ] + u[2, ]
    eps[1, ] <- Psi %*% eps0 + u[3, ]
    for (tt in 2:tmax) {
      eps[tt, ] <- Psi %*% eps[tt - 1, ] + u[tt + 2, ]
    }
    eps[(tmax - TT + 1) : tmax, ]
  })
  # C <- cov(t(sapply(epsilon, function(x) c(t(x)))))
  # sigma0 <- matrix(solve(diag(J^2) - kronecker(Psi, Psi), c(Sigma)), ncol = J)
  # C[1:3,1:3]
  # C[4:6,4:6]
  # C[7:9,7:9]
  # C[10:12,10:12]
  # C[13:15,13:15]
  errors <- do.call("rbind",epsilon)
  #lapply(1:nrow(epsilon), function(i)
  # matrix(epsilon[i, ], byrow = TRUE, ncol = J)))
  ytilde <- Xbeta + errors

  y.ord <- sapply(seq_len(J), function(j)
    cut(ytilde[, j], c(min(ytilde[, j]) - 1, c(theta[[j]]), max(ytilde[, j]) + 1),
        labels = FALSE), simplify = "array")

  df <- data.frame(firm_id  = rep(firm_id, J),
                   year_id  = rep(year_id, J),
                   outcome_id = rep(seq_len(J), each = I * TT),
                   response   = c(y.ord),
                   do.call("rbind", lapply(seq_len(J), function(i) X[, , j])))

  if (missings) {
    id_na <- sample(1:nrow(df), floor(nrow(df) * perc.miss))
    df[id_na, "rating"] <- NA
  }
  return(df)
}
set.seed(12345)

out <- list()
for (i in 1:100) {
  J  <- 3 ## number of responses
  TT <- 5 ## number of year
  I  <- 100 ## number of firms
  P  <- 2 ## number of covariates

  missings <- FALSE
  beta.matrix <- matrix(rep(c(2,-1),TT), byrow = TRUE, ncol = P)#cbind(c(2,1,3,-2),c(-1,1,0,1))#matrix(1, ncol = P, nrow = TT)#
  #matrix(rnorm(TT * P, 0, 2), ncol = P, byrow=TRUE)

  theta  <- list(c(- 1.0, 0, 1.5),
                 c( -0.5, 0, 0.5),
                 c(- 1.0, 0, 1))
  theta <- theta[1:J]


  Sigma <- diag(J)
  Sigma[upper.tri(Sigma)] <- Sigma[lower.tri(Sigma)] <-
    seq(0.95, 0.8, length.out = J * (J - 1)/2)


  #xi <-  0.9
  xiU <- c(0.8, 0.85, 0.9)
  xiL <- c(0.2, 0.25, 0.35)

  Psi <- xiU * diag(J)

  ## OR FULL PSI
  # repeat {
  # Psi <- matrix(rnorm(J^2), nrow = J)
  # if (all(abs(eigen(Psi)$values) < 1)) break
  # }

  data_toy_mvordflex2 <- simulate_3dim_1rep_for_loop(J = J, I = I, TT = TT, P = P,
                                                     beta.matrix = beta.matrix,
                                                     theta = theta,
                                                     Sigma = Sigma,
                                                     Psi = Psi, missings = missings,
                                                     perc.miss = perc.miss)
  res<- mvordflex(
    formula = MMO3(response, firm_id, year_id, outcome_id) ~ 0 + X1 + X2,
    data = data_toy_mvordflex2,
    error.structure = cor_MMO3(~1),
    coef.constraints = rep(1, q * TT),
    threshold.constraints = rep(1:q, TT),
    control = mvord::mvord.control(se = TRUE,
                                   solver = "newuoa",
                                   solver.optimx.control = list(maxit = 5000,
                                                                eval.max= 1000,
                                                                trace = 1)))
  out[[i]] <- summary(res)
}
save(out, file = "paper/simulation_I_100_T_5_id_1.rda")
