library("mvordflex")
data("data_toy_mvordflex", package = "mvordflex")
data("data_toy_mvordflex2", package = "mvordflex")
q <- 3 # number of multiple measurements
TT <- 5 # number of years

# Note that the number of responses is q*TT
## Probit
res<- mvordflex(
  formula = MMO3(response, firm_id, year_id, outcome_id) ~ 0 + X1 + X2,
  data = data_toy_mvordflex2,
  error.structure = cor_MMO3(~1),
  coef.constraints = rep(1, q * TT),
  threshold.constraints = rep(1:q, TT),
  control = mvord::mvord.control(se = FALSE,
                                 solver = "newuoa",
                                 solver.optimx.control = list(maxit = 5000,
                                                              eval.max= 1000,
                                                              trace = 1)))

summary(res)
Sigma
Psi

res<- mvordflex(
  formula = MMO3(response, firm_id, year_id, outcome_id) ~ 0 + X1 + X2,
  data = data_toy_mvordflex,
  error.structure = cor_MMO3(~1),
  coef.constraints = rep(1, q * TT),
  threshold.constraints = rep(1:q, TT), PL.lag = 1)


## Logit
res_logit <- mvordflex(
  formula = MMO3(response, firm_id, year_id, outcome_id) ~ 0 + X1 + X2,
  data = data_toy_mvordflex,
  link = mvord::mvlogit(),
  error.structure = cor_MMO3(~1),
  coef.constraints = rep(1, q * TT),
  threshold.constraints = rep(1:q, TT))


print(res)
summary(res)
mvord::thresholds(res)
coefficients(res)
head(error_structure(res))

# Note that the number of responses is q*TT
res_ar1 <- mvordflex(
  formula = MMO3(response, firm_id, year_id, outcome_id) ~ 0 + X1 + X2,
  data = data_toy_mvordflex,
  error.structure = cor_MMO3_ar1(~1),
  coef.constraints = rep(1, q * TT),
  threshold.constraints = rep(1:q, TT), PL.lag = 1)

# Note that the number of responses is q*TT
res_cross <- mvordflex(
  formula = MMO3(response, firm_id, year_id, outcome_id) ~ 0 + X1 + X2,
  data = data_toy_mvordflex,
  error.structure = cor_MMO3_cross(~1),
  coef.constraints = rep(1, q * TT),
  threshold.constraints = rep(1:q, TT), PL.lag = 1)

AIC(res_ar1, res_cross, res)

## Probit
res <- mvordflex(
  formula = MMO3(response, firm_id, year_id, outcome_id) ~ 0 + X1 + X2,
  data = data_toy_mvordflex2,
  error.structure = cor_MMO3(~1, Psi.diag = FALSE),
  #coef.constraints = rep(1, q * TT),
  #threshold.constraints = rep(1:q, TT),
  PL.lag = 1,
  control = mvord::mvord.control(se = FALSE))
,
                                 solver = "newuoa",
                                 solver.optimx.control = list(maxit = 5000,
                                                              eval.max= 1000,
                                                              trace = 1)))

summary(res)
res$rho$optRes
