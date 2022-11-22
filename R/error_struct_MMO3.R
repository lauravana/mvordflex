#' @title Error Structures in mvordflex
#' @description Different \code{error.structures} are available in \pkg{mvordflex}.
#' For three dimensional data in particular we have
#' \itemize{
#' \item correlation structure based on a multivariate AR(1) process \code{cor_MMO3(~ 1)}
#' \item correlation structure for three dimensional panel data where only
#'       cross-sectional correlation is accounted for \code{cor_MMO3_cross(~ 1)}
#' \item correlation structure for three dimensional panel data where only longitudinal
#'       correlation is accounted for by an AR(1) structure \code{cor_MMO3_ar1(~ 1)}
#' }
#' @param formula \code{\link{formula}} object
#' @param value specifies values of the correlation (and variance) parameters.
#' specifies values of the correlation (and variance) parameters. For \code{cor_MMO3()}
#' the first parameters correspond to the cross-sectional error structure while the
#' latter correspond to the AR(1) coefficients for each outcome. Default is \code{value
#' = numeric(0)} object. In this case the correlation parameters are initialized with zero.
#' @param fixed logical specifying whether the parameters of the error structure should not be optimized in the procedure, but will be \cr
#' fixed to the values specified in the argument \code{value}. Defaults to \code{fixed = FALSE}.
#' @export
#' @name error_struct
#' @export
cor_MMO3 <-
  ## Constructor for the cor_MMO3 class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE, Psi.diag = TRUE)
  {
    obj <- list(name = "cor_MMO3",
                formula = formula,
                type = "correlation", value = value, fixed = fixed,
                Psi.diag = Psi.diag)
    attr(obj, "formula") <- formula
    class(obj) <- c("cor_MMO3", "error_struct")
    obj
  }

#' @rdname error_struct
#' @export
cor_MMO3_cross <-
  ## Constructor for the cor_MMO3_cross class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE, Psi.diag = TRUE)
  {
    obj <- list(name = "cor_MMO3_cross",
                formula = formula,
                type = "correlation", value = value, fixed = fixed,
                Psi.diag = Psi.diag)
    attr(obj, "formula") <- formula
    class(obj) <- c("cor_MMO3_cross", "error_struct")
    obj
  }

#' @rdname error_struct
#' @export
cor_MMO3_ar1 <-
  ## Constructor for the cor_MMO3_ar1 class
  function(formula = ~ 1, value = numeric(0), fixed = FALSE, Psi.diag = TRUE)
  {
    obj <- list(name = "cor_MMO3_ar1",
                formula = formula,
                type = "correlation", value = value, fixed = fixed,
                Psi.diag = Psi.diag)
    attr(obj, "formula") <- formula
    class(obj) <- c("cor_MMO3_ar1", "error_struct")
    obj
  }
########################################

build_error_struct <-
  ## extractor for correlation matrix
  function(eobj, ...) UseMethod("build_error_struct")

build_error_struct_fixed <-
  ## extractor for correlation matrix
  function(eobj, ...) UseMethod("build_error_struct_fixed")

start_values <-
  ## extractor for correlation matrix
  function(eobj, ...) UseMethod("start_values")

initialize <-
  ## initializes the structures
  function(eobj, ...) UseMethod("initialize")

finalize_fun <-
  ## finalizes the structures
  function(eobj, ...) UseMethod("finalize_fun")

finalize <-
  ## finalizes the structures
  function(eobj, ...) UseMethod("finalize")

get_covariate <-
  ## initializes the structures
  function(eobj, ...) UseMethod("get_covariate")

init_fun <-
  ## initializes the structures
  function(eobj, ...) UseMethod("init_fun")

#################
##   Methods for error_struct
#################
formula.error_struct <-
  ## Accessor for the covariate formula
  function(x, ...) eval(attr(x, "formula"))

#data.x = data$x
get_covariate.error_struct <- function(eobj, data.x, contrasts) {
  covar_mat <- lapply(data.x, function(x)
    suppressWarnings(model.matrix(formula(eobj),
                                  model.frame(formula(eobj), x, na.action = function(x) x),
                                  contrasts.arg = contrasts)))
  ## check if covariate matrices are equal
  if (!all(sapply(1:(length(covar_mat) - 1), function(i)
    all(covar_mat[[i]] == covar_mat[[i+1]], na.rm = T)))) {
    stop("Covariates in error structure must be
         constant across outcomes!")
  }
  # make one matrix
  covar_mat1 <- sapply(1:ncol(covar_mat[[1]]), function(k){
    xtcol <- do.call(cbind,lapply(covar_mat, `[`,, k))
    xtcol_final <- apply(xtcol,1,function(i) unique(i[!is.na(i)]))
    xtcol_final
  })
  attributes(covar_mat1) <- attributes(covar_mat[[1]])[1:2]
  covar_mat1
}

initialize.error_struct <-
  ## initializes some attributes of error_struct objects
  ## takes as data the output on mvord_data
  function(eobj, data, contrasts)
  {
    attr(eobj, "ynames") <- colnames(data$y)
    #attr(eobj, "subjnames") <- rownames(data$y)
    attr(eobj, "ndim") <- length(data$x)
    attr(eobj, "nobs") <- nrow(data$y)
    attr(eobj, "covariate") <-
      get_covariate(eobj, data.x = data$x, contrasts = contrasts)
    eobj
  }

update.error_struct <-
  ## initializes some attributes of error_struct objects
  ## takes as data the output on mvord_data
  function(eobj, data, contrasts)
  {
    attr(eobj, "ynames") <- colnames(data$y)
    attr(eobj, "ndim") <- length(data$x)
    attr(eobj, "nobs") <- nrow(data$y)
    attr(eobj, "covariate") <-
      get_covariate(eobj, data.x = data$x, contrasts = contrasts)
    eobj
  }

finalize.error_struct <-
  ## initializes some attributes of error_struct objects
  ## takes as data the output on mvord_data
  function(eobj, tpar)
  {
    eobj <- finalize_fun(eobj, tpar)
    eobj$value_tmp <- NULL
    attr(eobj, "npar.cor") <- NULL
    attr(eobj, "npar.sd") <- NULL
    eobj
  }


initialize.MMO3 <-
  ## initializes some attributes of error_struct objects
  ## takes as data the output on mvord_data
  function(eobj, data, contrasts)
  {
    attr(eobj, "ynames") <- colnames(data$y)
    attr(eobj, "ndim") <- length(data$x)
    attr(eobj, "nobs") <- nrow(data$y)
    attr(eobj, "covariate") <-
      get_covariate(eobj, data.x = data$x, contrasts = contrasts)
    eobj
  }

### Function for creating the (JxTT) x (JxTT) correlation matrix
# #' @param blocks a list of matrices
# #'              Sigma, Psi * Sigma, Psi^2 * Sigma, ...,Psi^(T-1) * Sigma
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

make_Q_cor_MMO3 <- function(tpar, eobj){
  ndim <- attr(eobj, "ndim")
  ndim_j <- attr(eobj, "ndim_j")
  ndim_t <- attr(eobj, "ndim_t")
  npar_sigma <- attr(eobj, "npar_sigma")
  npar_psi <- attr(eobj, "npar_psi")
  covar <- attr(eobj, "covariate")
  tpar_sigma <- tpar[seq_len(npar_sigma)]
  tpar_psi   <- tpar[seq_len(npar_psi) + npar_sigma]

  ## Sigma
  angles <- pi * exp(tpar_sigma)/(1 + exp(tpar_sigma))
  cosmat <- diag(ndim_j)
  cosmat[lower.tri(cosmat)] <- if(length(angles)>0) cos(angles) else 0
  S1 <- matrix(0, nrow = ndim_j, ncol = ndim_j)
  S1[, 1L] <- 1
  S1[lower.tri(S1, diag = T)][-(1:ndim_j)] <- if(length(angles)>0) sin(angles) else 1
  #S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
  tLmat <- sapply(seq_len(ndim_j),
                  function(j) cosmat[j, ] * cumprod(S1[j, ]))
  sigma <- crossprod(tLmat)


  # PSI matrix
  if (npar_psi == 0) {
    psi <- diag(z2r(rep(0, ndim_j)))
    ar_blocks <- lapply(0:(ndim_t - 1), function(t) {
      s_ar <- if (t == 0) diag(ndim_j) else psi ^ t
      crossprod(s_ar, sigma)
    })
  }
  if(eobj$Psi.diag & npar_psi != 0) {
    psi <- diag(z2r(tpar_psi))
    ar_blocks <- lapply(0:(ndim_t - 1), function(t) {
      s_ar <- if (t == 0) diag(ndim_j) else psi ^ t
      crossprod(s_ar, sigma)
    })
  }
  if(!eobj$Psi.diag & npar_psi != 0) {
    psi <-  make_stationary_psi(tpar_psi)
    ar_blocks <- vector("list", length = ndim_t)
    for (tt in 0:(ndim_t - 1)) {
      s_ar <- if (tt == 0) diag(ndim_j) else s_ar %*% psi
      ar_blocks[[tt + 1]] <- s_ar %*% sigma
    }
  }

  S <- toeplitz.block(ar_blocks)
  S
}

error_structure.cor_MMO3 <- function(eobj, type, ...){
  par <- attr(eobj, "par")
  npar <- attr(eobj, "npar")# npar <- length(par)
  ndim <- attr(eobj, "ndim")
  ndim_t <- attr(eobj, "ndim_t")
  ndim_j <- attr(eobj, "ndim_j")
  covar <- attr(eobj, "covariate")
  ynames <- attr(eobj, "ynames")
  nlev <- NCOL(covar)
  npar.cor <- npar/nlev
  npar_sigma <- attr(eobj, "npar_sigma")
  npar_psi <- attr(eobj, "npar_psi")

  par_sigma <- par[seq_len(npar_sigma)]
  par_psi <- par[seq_len(npar_psi)+ npar_sigma]


  S1 <- diag(ndim_j)
  S1[lower.tri(S1)] <- par_sigma[(l - 1) * npar1 + seq_len(npar1)]
  S1[upper.tri(S1)] <- t(S1)[upper.tri(S1)]
  sigma <- S1

  if (eobj$Psi.diag){
    psi <- diag(z2r(par_psi))
    ar_blocks <- lapply(0:(ndim_t - 1), function(t) {
      s_ar <- if (t == 0) diag(ndim_j) else psi ^ t
      crossprod(s_ar, sigma)
    })
  } else {
    psi <-  make_stationary_psi(par_psi)
    ar_blocks <- vector("list", length = ndim_t)
    for (tt in 0:(ndim_t - 1)) {
      s_ar <- if (tt == 0) diag(ndim_j) else s_ar %*% psi
      ar_blocks[[tt + 1]] <- s_ar %*% sigma
    }
  }
  S <- toeplitz.block(ar_blocks)
  S
}

error_structure.cor_MMO3_cross <- function(eobj, type, ...){
  par <- attr(eobj, "par")
  npar <- attr(eobj, "npar")# npar <- length(par)
  ndim <- attr(eobj, "ndim")
  ndim_t <- attr(eobj, "ndim_t")
  ndim_j <- attr(eobj, "ndim_j")
  covar <- attr(eobj, "covariate")
  ynames <- attr(eobj, "ynames")
  nlev <- NCOL(covar)
  npar.cor <- npar/nlev
  npar_sigma <- attr(eobj, "npar_sigma")
  npar_psi <- 0# attr(eobj, "npar_psi")

  par_sigma <- par[seq_len(npar_sigma)]
  par_psi <- par[seq_len(npar_psi)+ npar_sigma]


  S1 <- diag(ndim_j)
  S1[lower.tri(S1)] <- par_sigma
  S1[upper.tri(S1)] <- t(S1)[upper.tri(S1)]
  sigma <- S1


  ar_blocks <- sigma

  S <- toeplitz.block(ar_blocks)
  S
}

error_structure.cor_MMO3_ar1 <- function(eobj, type, ...){
  par <- attr(eobj, "par")
  npar <- attr(eobj, "npar")# npar <- length(par)
  ndim <- attr(eobj, "ndim")
  ndim_t <- attr(eobj, "ndim_t")
  ndim_j <- attr(eobj, "ndim_j")
  covar <- attr(eobj, "covariate")
  ynames <- attr(eobj, "ynames")
  nlev <- NCOL(covar)
  npar.cor <- npar/nlev
  npar_sigma <- attr(eobj, "npar_sigma")
  npar_psi <- attr(eobj, "npar_psi")
  par_sigma <- par[seq_len(npar_sigma)]
  par_psi   <- par[seq_len(npar_psi) + npar_sigma]
  if (eobj$Psi.diag){
    psi <- diag(z2r(par_psi))
    ar_blocks <- lapply(0:(ndim_t - 1), function(t) {
      s_ar <- if (t == 0) diag(ndim_j) else psi ^ t
      crossprod(s_ar, diag(ndim_j))
    })
  } else {
    psi <- make_stationary_psi(par_psi)
    ar_blocks <- vector("list", length = ndim_t)
    for (tt in 0:(ndim_t - 1)) {
      s_ar <- if (tt == 0) diag(ndim_j) else s_ar %*% psi
      ar_blocks[[tt + 1]] <- s_ar %*% sigma
    }
  }
  S <- toeplitz.block(ar_blocks)
  S
}


##############################
#### Methods for cor_MMO3 ####
##############################
# eobj <- rho$error.structure
# y <- rho$y
start_values.cor_MMO3 <- start_values.cor_MMO3_cross <- start_values.cor_MMO3_ar1 <-
  function(eobj, y = NULL) {
    ## builds starting values for the cor_MMO3 structures
    rep(0, attr(eobj, "npar"))
}


init_fun.cor_MMO3 <-
  function(eobj,  data, contrasts)
  {
    ## initializes some attributes of cor_general eobjs
    form <- formula(eobj)
    if (length(all.vars(form)) > 1)
      stop("Only one factor is supported in cor_general.")
    ## if intercept included rewrite formula without
    if (length(all.vars(form)) == 1 & attr(terms(form), "intercept") == 1)
      attr(eobj, "formula") <- as.formula(sprintf("~ 0 + %s", all.vars(form)))
    eobj <- initialize.MMO3(eobj, data, contrasts)

    n    <- attr(eobj, "nobs")
    ndim <- attr(eobj, "ndim")
    ndim_t <- attr(eobj, "ndim_t")
    ndim_j <- attr(eobj, "ndim_j")

    npar_sigma <-  ndim_j * (ndim_j - 1)/2
    npar_psi <-  ifelse(eobj$Psi.diag, ndim_j, ndim_j ^ 2)

    attr(eobj, "npar_sigma") <- npar_sigma
    attr(eobj, "npar_psi") <- npar_psi
    npar <- npar_sigma + npar_psi
    attr(eobj, "npar") <- npar

    r <- eobj$value

    if (is.null(eobj$fixed)) eobj$fixed  <- FALSE
    attr(eobj, "npar.cor") <- ifelse(eobj$fixed, 0, attr(eobj, "npar"))
    attr(eobj, "npar.sd") <- 0
    attr(eobj, "npar") <-   attr(eobj, "npar.cor") + attr(eobj, "npar.sd")
    if(length(all.vars(form)) == 1 && !is.factor(data$x[[1]][, all.vars(form)]))
      stop("For cor_MMO3 covariate must be factor!")
    eobj
  }

init_fun.cor_MMO3_cross <-
  function(eobj,  data, contrasts)
  {
    ## initializes some attributes of cor_general eobjs
    form <- formula(eobj)
    if (length(all.vars(form)) > 1)
      stop("Only one factor is supported in cor_general.")
    ## if intercept included rewrite formula without
    if (length(all.vars(form)) == 1 & attr(terms(form), "intercept") == 1)
      attr(eobj, "formula") <- as.formula(sprintf("~ 0 + %s", all.vars(form)))
    eobj <- initialize.MMO3(eobj, data, contrasts)

    n    <- attr(eobj, "nobs")
    ndim <- attr(eobj, "ndim")
    ndim_t <- attr(eobj, "ndim_t")
    ndim_j <- attr(eobj, "ndim_j")

    npar_sigma <- ndim_j * (ndim_j - 1)/2
    npar_psi   <- 0
    attr(eobj, "npar_sigma") <- npar_sigma
    attr(eobj, "npar_psi") <- npar_psi
    npar <- npar_sigma + npar_psi
    attr(eobj, "npar") <- npar

    r <- eobj$value
    if (is.null(eobj$fixed)) eobj$fixed  <- FALSE
    attr(eobj, "npar.cor") <- ifelse(eobj$fixed, 0,   attr(eobj, "npar") )
    attr(eobj, "npar.sd") <- 0
    attr(eobj, "npar") <-   attr(eobj, "npar.cor") + attr(eobj, "npar.sd")
    if(length(all.vars(form)) == 1 && !is.factor(data$x[[1]][, all.vars(form)]))
      stop("For cor_MMO3_cross covariate must be factor!")
    eobj
  }

init_fun.cor_MMO3_ar1 <-
  function(eobj,  data, contrasts)
  {
    ## initializes some attributes of cor_general eobjs
    form <- formula(eobj)
    ## if intercept included rewrite formula without
    if (length(all.vars(form)) == 1 & attr(terms(form), "intercept") == 1)
      attr(eobj, "formula") <- as.formula(sprintf("~ 0 + %s", all.vars(form)))
    eobj <- initialize.MMO3(eobj, data, contrasts)

    n    <- attr(eobj, "nobs")
    ndim <- attr(eobj, "ndim")
    ndim_t <- attr(eobj, "ndim_t")
    ndim_j <- attr(eobj, "ndim_j")

    npar_sigma <- 0
    npar_psi   <- ifelse(eobj$Psi.diag, ndim_j, ndim_j ^ 2)
    attr(eobj, "npar_sigma") <- npar_sigma
    attr(eobj, "npar_psi") <- npar_psi
    npar <- npar_sigma + npar_psi
    attr(eobj, "npar") <- npar

    r <- eobj$value
    # ## check value
    # if (length(r) == 0) r <- numeric(ndim *  (ndim - 1)/2) # if default set to zero
    # if (length(r) == n) r <- cbind(r)
    # if (length(r) == npar_psi *  (npar_psi - 1)/2)  r <- matrix(rep.int(r, n), ncol = npar_psi * (npar_psi - 1)/2, byrow= T) # if only one vector of length npar_psi *  (npar_psi - 1)/2 default set to zero
    #
    # if (nrow(r) != n) stop("Number of rows of argument value in cor_general() is not equal to number of subjects.")
    #
    # ## TODO - check positive semi-definiteness??
    # ## end check
    # eobj$value_tmp <- r
    if (is.null(eobj$fixed)) eobj$fixed  <- FALSE
    attr(eobj, "npar.cor") <- ifelse(eobj$fixed, 0, npar)
    attr(eobj, "npar.sd") <- 0
    attr(eobj, "npar") <-   attr(eobj, "npar.cor") + attr(eobj, "npar.sd")
    if(length(all.vars(form)) == 1 && !is.factor(data$x[[1]][, all.vars(form)]))
      stop("For cor_MMO3_ar1 covariate must be factor!")
    eobj
  }

# build_error_struct_fixed.cor_MMO3 <-
#   ## builds the correlation matrix when fixed = T of cor_general objects
#   function(eobj, tpar = NULL)
#   {
#     ## takes the transformed parameters and builds initializes some attributes of cor_general objects
#     sd <- rep.int(1, attr(eobj, "ndim"))
#     return(list(rVec = eobj$value_tmp, sdVec = sd))
#   }


# eobj <- rho$error.structure
# tpar <- par_sigma
#tpar <- c(1,2,3)

build_error_struct.cor_MMO3 <- build_error_struct.cor_MMO3_cross <-
  build_error_struct.cor_MMO3_ar1 <-
  function(eobj, tpar)
  {
    ## takes the transformed parameters and builds initializes some attributes of cor_general eobjs
    ndim <- attr(eobj, "ndim")
    covar <- attr(eobj, "covariate")

    Q <- make_Q_cor_MMO3(tpar, eobj)

    corr_pars <- Q[lower.tri(Q)]#sapply(Q, function(x) x[lower.tri(x)])

    rVec <- tcrossprod(covar, corr_pars)#tcrossprod(covar, s[lower.tri(s)])#tcrossprod(covar, corr_pars)
    sd <- rep(1, ndim)
    return(list(rVec = rVec, sdVec = sd))
  }


finalize_fun.cor_MMO3 <- finalize_fun.cor_MMO3_cross <- finalize_fun.cor_MMO3_ar1 <- function(eobj, tpar) {
  if (eobj$fixed){
    attr(eobj, "par") <- numeric(0)
  } else {
    ndim <- attr(eobj, "ndim")
    ndim_j <- attr(eobj, "ndim_j")
    ndim_t <- attr(eobj, "ndim_t")
    npar_sigma <- attr(eobj, "npar_sigma")
    npar_psi <- attr(eobj, "npar_psi")
    covar <- attr(eobj, "covariate")
    tpar_sigma <- tpar[seq_len(npar_sigma)]
    tpar_psi   <- tpar[seq_len(npar_psi)+ npar_sigma]

    angles <- pi * exp(tpar_sigma)/(1 + exp(tpar_sigma))
    cosmat <- diag(ndim_j)
    cosmat[lower.tri(cosmat)] <- if(length(angles)>0) cos(angles) else 0
    S1 <- matrix(0, nrow = ndim_j, ncol = ndim_j)
    S1[, 1L] <- 1
    S1[lower.tri(S1, diag = T)][-(1:ndim_j)] <- if(length(angles)>0) sin(angles) else 1
    #S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
    tLmat <- sapply(seq_len(ndim_j),
                    function(j) cosmat[j, ] * cumprod(S1[j, ]))
    sigma <- crossprod(tLmat)


    if (eobj$Psi.diag) {
      tpar_psi1 <- if (npar_psi==0) rep(0, ndim_j) else tpar_psi
      psi <- diag(z2r(tpar_psi1))
    } else {
      tpar_psi1 <- if (npar_psi==0) rep(0, ndim_j^2) else tpar_psi
      psi <- make_stationary_psi(tpar_psi1)
    }
    ## names
    ynames <- attr(eobj, "ynames")
    ind <- combn(ndim_j, 2)
    if (eobj$formula == ~1) {
      ## correlation names
      names.sigma <-
        apply( combn(ndim_j, 2), 2,  function(x)
          sprintf("corr %s %s", x[1], x[2]))
      names.psi <- if (eobj$Psi.diag) paste0("psi ", 1:ndim_j) else c(outer(1:m, 1:m, function(x,y) paste("psi", x, y)))

    }

    par_sigma <- if (npar_sigma != 0) sigma[lower.tri(sigma)] else NULL
    par_psi <- NULL
    if (npar_psi != 0 & eobj$Psi.diag)  par_psi <- diag(psi)
    if (npar_psi != 0 & !eobj$Psi.diag) par_psi <- c(psi)

    attr(eobj, "par") <- c(par_sigma, par_psi)
    attr(eobj, "parnames") <-   c(if (npar_sigma != 0) names.sigma,
                                  if (npar_psi   != 0) names.psi)
  }
  eobj
}
sqrootmat  <- function(M) {
  e <- eigen(M)
  e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
}

make_stationary_psi <- function(tpar) {
  ## Works for now only for symmetric Psi
  m <- sqrt(length(tpar)) # this is ndim_j
  l <- tpar[seq_len(m * (m - 1)/2)]                     ## cholesky of V positive definite
  d <- tpar[m * (m - 1)/2 + seq_len(m)]                 ## diag V positive definite
  s <- tpar[m * (m - 1)/2 + m + seq_len(m * (m - 1)/2)] ## Q orthogonal
  ## Transformations
  ### For V - positive definite: LDL'
  logV <- matrix(0, nrow = m, ncol = m)
  logV[upper.tri(logV)] <- l
  logV <- logV + t(logV)
  diag(logV) <- d
  e <- eigen(logV)
  U <- e$vectors
  loglambda <- e$values
  V <- U %*% diag(exp(loglambda)) %*% t(U)
  ### For Q - orthogonal - for now not working
  S <- matrix(0, nrow = m, ncol = m) # skew symmetric S' = -S
  #S <- diag(m)
  S[lower.tri(S)] <- s
  S <- S - t(S)
  Im <- diag(m)
  ## Givens parametrization
  # index <- 1
  # Q <- diag(m)
  # for(i in 1:(m-1)){
  #   for(j in (i+1):m){
  #     cat(i,j)
  #     if (i == 1) {
  #       theta <- 2 * pi * (1 - exp(s[index]))/(1 + exp(s[index]))
  #     } else {
  #       theta <- pi * (1 - exp(s[index]))/(1 + exp(s[index]))
  #     }
  #
  #     cos1 <- cos(theta);
  #     sin1 <- sin(theta);
  #     Qi_tmp <- cos1 * Q[,i] - sin1 * Q[,j];
  #     Qj_tmp <- sin1 * Q[,i] + cos1 * Q[,j];
  #     Q[, i] <- Qi_tmp;
  #     Q[, j] <- Qj_tmp;
  #     index <- index + 1;
  #   }
  #   #index <- 1;
  # }
  ind <- combn(m, 2)
  Q <- diag(m)
  for (k in 1:ncol(ind)){
    i <- ind[1, k]
    j <- ind[2, k]
    theta <-  ifelse(i == 1, 2, 1) * pi * (1 - exp(s[k]))/(1 + exp(s[k]))
    Gk <- diag(m)
    Gk[i, i] <- Gk[j, j] <- cos(theta)
    Gk[i, j] <-   sin(theta)
    Gk[j, i] <- - sin(theta)
    Q <- Q %*% Gk
  }
  #Q <- PL.Q(Im, S)
  #Q <-  diag(c(1, 1,-1)) %*% (Im - S) %*% solve(Im + S)
  # Q <- solve(Im - S, Im + S)
  #Q <- Im#(Im - S) %*% solve(Im+S)
  ### make matrix A = V^{1/2}Q(V + M)^{-1/2}, M-fixed
  Vhalf <- U %*% diag(exp(loglambda/2)) %*% t(U) #
  psi   <- Vhalf %*% Q %*% sqrootmat(chol2inv(chol(V + Im)))
  psi
}
