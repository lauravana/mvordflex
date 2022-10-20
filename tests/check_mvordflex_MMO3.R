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

simulate_3dim_1rep <- function(J, I, TT, P,
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
  ## Correlation structure for Psi = diag(xi_j)
  ar_blocks <- lapply(0:(TT - 1), function(t) {
    s_ar <- if (t == 0) diag(J) else Psi ^ t
    crossprod(s_ar, Sigma)
  })
  S <- toeplitz.block(ar_blocks)
  ## Simulate errors
  epsilon <- mvnfast::rmvn(I, rep(0, J * TT), S)
  errors <- do.call("rbind", lapply(1:nrow(epsilon), function(i)
    matrix(epsilon[i, ], byrow = TRUE, ncol = J)))

  ytilde <- Xbeta + errors

  y.ord <- sapply(seq_len(J), function(j)
    cut(ytilde[, j], c(min(ytilde[, j]) - 1, c(theta[[j]]), max(ytilde[, j]) + 1),
        labels = FALSE), simplify = "array")

  df <- data.frame(firm_id  = rep(firm_id, J),
                   year_id  = rep(year_id, J),
                   rater_id = rep(seq_len(J), each = I * TT),
                   rating   = c(y.ord),
                   do.call("rbind", lapply(seq_len(J), function(i) X[, , j])))

  if (missings) {
    id_na <- sample(1:nrow(df), floor(nrow(df) * perc.miss))
    df[id_na, "rating"] <- NA
  }
  return(df)
}

### Simulate and estimate mvordflex nrep times
#' @param J number of ordinal responses
#' @param I number of subjects or firms
#' @param TT number of time points
#' @param P number of covariates to be generated from a standard normal (we assume them constant over the responses)
#' @param beta.matrix a TT x P matrix of coefficients which are allowed to vary over time
#' @param theta a list of length J where each element contains a vector of threshold parameters for response j
#' @param Sigma a J x J matrix which corresponds to the lag 0 correlation among the J responses
#' @param Psi a J x J  matrix which in the simplest case is diagonal containing the autoregressive coefficients for each rater.
#'        Otherwise it could be a matrix corresponding to a VAR type of model.
#' @param missings logical indicating whether missing values should be generated at random
#' @param perc.miss percentage of missing values which is constant over all responses

simulate_3dim <- function(nrep, J, I, TT, P,
                          beta.matrix = NULL, theta = NULL, Sigma = NULL, Psi = NULL,
                          missings = FALSE, perc.miss = 0.2) {
  sapply(seq_len(nrep), function(i) {
    df <- simulate_3dim_1rep(J = J, I = I, TT = TT, P = P,
                             beta.matrix = beta.matrix,
                             theta = theta,
                             Sigma = Sigma,
                             Psi = Psi, missings = missings,
                             perc.miss = perc.miss)
    res <- mvord(formula = MMO3(rating, firm_id, year_id, rater_id) ~ 0 + X1 + X2,
                 data = df, error.structure = cor_MMO3(~1),
                 coef.constraints = rep(1:TT, each = J),
                 threshold.constraints = rep(1:J, TT))
    s <- summary(res)
    return(s)

  })

}



## Design of simulation study
## 1000 rep
## n=1000
## J = 3, T = 5

n_vec <- c(1000)

set.seed(12345)

J  <- 3 ## number of responses
TT <- 4 ## number of year
I  <- 100 ## number of firms
P  <- 2 ## number of covariates

missings <- FALSE
beta.matrix <- #cbind(c(2,1,3,-2),c(-1,1,0,1))#matrix(1, ncol = P, nrow = TT)#
  matrix(rnorm(TT * P, 0, 2), ncol = P, byrow=TRUE)

theta  <- list(c(- 1.0, 0, 1.5),
               c( -0.5, 0, 0.5),
               c(- 1.0, 0, 1))
theta <- theta[1:J]


Sigma <- diag(J)
Sigma[upper.tri(Sigma)] <- Sigma[lower.tri(Sigma)] <-
  seq(0.95, 0.8, length.out = J * (J - 1)/2)


#xi <-  0.9
xi <- round(runif(J, 0.8, 0.9), 2)
Psi <- xi * diag(J)


df_3dim <- simulate_3dim_1rep(J = J, I = 1000, TT = TT, P = P,
                              beta.matrix = beta.matrix, theta = theta,
                              Sigma = Sigma, Psi = Psi, missings = FALSE, perc.miss = 0)


library(mvordflex)
res_MMO3_test <- mvord(formula = MMO3(rating, firm_id, year_id, rater_id) ~ 0 + X1 + X2,
                       data = (df_3dim),
                       error.structure = cor_MMO3(~1),
                       coef.constraints = rep(1:TT, each = J),
                       threshold.constraints = rep(1:J,TT))
res.summary <- summary(res_MMO3_test)
cbind(Estimate = c(res.summary$thresholds$Estimate,
                   res.summary$coefficients$Estimate,
                   res.summary$error.structure$Estimate),
      True = c(unlist(theta), c(beta.matrix), Sigma[upper.tri(Sigma)], xi))

rho <- res_MMO3_test$rho


# paste(format(res.summary$thresholds$Estimate), collapse = ",")
# paste(format(res.summary$coefficients$Estimate), collapse = ",")
# paste(format(res.summary$error.structure$Estimate), collapse = ",")
mvord:::check(all.equal(res.summary$thresholds$Estimate, c(-0.988159740166693, 0.00438259573864504, 1.55021869952464,
                                                           -0.519125525471071, -0.0120591404810041, 0.500569066851932,
                                                           -0.987847466599132, 0.0202335070082556, 1.02103588820272), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$Estimate, c(1.21068036395154, -0.208798249535033, 1.17349676315999, 1.2218089317592,
                                                             1.47846684950275, -0.943989437508043, -3.61538830543199, -0.513253507216111), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$Estimate, c(0.955035462031249, 0.880720689267029, 0.807421877607695, 0.831848043092694, 0.828977428641756, 0.815718915373051), tolerance = tolerance))
mvord:::check(all.equal(res.summary$thresholds$`Std. Error`, c(0.0399189784939852, 0.0350222836494116, 0.0452354237741449,
                                                               0.0368978788200778, 0.0345642181745102, 0.0354708824440899,
                                                               0.0405392865593343, 0.0348651573685353, 0.0385077589394305), tolerance = tolerance))
mvord:::check(all.equal(res.summary$coefficients$`Std. Error`, c(0.0433975823345028, 0.0259795048857014, 0.0485461781375225, 0.0405537988459686,
                                                                 0.0490195753399894, 0.0332631024042478, 0.111395365065635, 0.0300579052219538), tolerance = tolerance))
mvord:::check(all.equal(res.summary$error.structure$`Std. Error`, c(0.00429683238466654, 0.00874602660742379, 0.013378771719646, 0.0111508616105756, 0.0131829835861842, 0.0131168355262699), tolerance = tolerance))
mvord:::check(all.equal(logLik(res_MMO3_test)[[1]], -97447.22, tolerance = tolerance))
mvord:::check(all.equal(AIC(res_MMO3_test), 195503.9, tolerance = tolerance))
mvord:::check(all.equal(BIC(res_MMO3_test), 196999.5, tolerance = tolerance))


## Check PL.lag
library(mvordflex)
res_MMO3_PL <- mvord(formula = MMO3(rating, firm_id, year_id, rater_id) ~ 0 + X1 + X2,
                       data = (df_3dim),
                       error.structure = cor_MMO3(~1),
                       coef.constraints = rep(1:TT, each = J),
                       threshold.constraints = rep(1:J,TT),
                     PL.lag = 2L)
res.summary <- summary(res_MMO3_PL)
cbind(Estimate = c(res.summary$thresholds$Estimate,
                   res.summary$coefficients$Estimate,
                   res.summary$error.structure$Estimate),
      True = c(unlist(theta), c(beta.matrix), Sigma[upper.tri(Sigma)], xi))
