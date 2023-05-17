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
  function(formula = ~ 1, value = numeric(0), fixed = FALSE)
  {
    obj <- list(name = "cor_MMO3_cross",
                formula = formula,
                type = "correlation", value = value, fixed = fixed)
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
    #attr(eobj, "subjnames") <- rownames(data$y)
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

    #  attr(eobj, "subjnames") <- NULL
    #attr(eobj, "ynames") <- NULL
    #  attr(eobj, "ndim") <- NULL
    # attr(eobj, "nobs") <- NULL
    # attr(eobj, "covariate") <- NULL
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
    #attr(eobj, "subjnames") <- rownames(data$y)
    attr(eobj, "ndim") <- length(data$x)
    # attr(eobj, "npar_sigma") <- length(data$x)
    # attr(eobj, "npar_psi") <- length(data$x)
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
  nlev <- 1
  npar1 <- attr(eobj, "npar_sigma")
  npar1_psi <- attr(eobj, "npar_psi")
  tpar_sigma <- tpar[seq_len(npar_sigma)]
  tpar_psi <- tpar[seq_len(npar_psi)+ npar_sigma]

  nu <- tpar_sigma
  angles <- pi * exp(nu)/(1 + exp(nu))
  cosmat <- diag(ndim_j)
  cosmat[lower.tri(cosmat)] <- if(length(angles)>0) cos(angles) else 0
  S1 <- matrix(0, nrow = ndim_j, ncol = ndim_j)
  S1[, 1L] <- 1
  S1[lower.tri(S1, diag = T)][-(1:ndim_j)] <- if(length(angles)>0) sin(angles) else 1
  #S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
  tLmat <- sapply(seq_len(ndim_j),
                    function(j) cosmat[j, ] * cumprod(S1[j, ]))
  sigma <- crossprod(tLmat)

  tpar_psi1 <- if(npar1_psi == 0) rep(0, ndim_j) else tpar_psi
  psi <- diag(z2r(tpar_psi1))

  # sigma0 <-
  #  matrix(solve(diag(ndim_j^2) - kronecker(psi, psi),
  #               c(sigma)),
  #         ncol = ndim_j)
  sigma0_vec <- 1/(1 - diag(kronecker(psi, psi))) * c(sigma)
  sigma0 <- matrix(sigma0_vec, ncol = ndim_j)

  ar_blocks <- lapply(0:(ndim_t - 1), function(t) {
    s_ar <- if (t == 0) diag(ndim_j) else psi ^ t
    tcrossprod(sigma0, s_ar)
  })

  S <- list(toeplitz.block(ar_blocks))
  S
}

error_structure.cor_MMO3 <- function(eobj, type, ...){
  par <- attr(eobj, "par")
  npar <- attr(eobj, "npar")
  ndim_t <- attr(eobj, "ndim_t")
  ndim_j <- attr(eobj, "ndim_j")
  npar_sigma <- attr(eobj, "npar_sigma")
  npar_psi <- attr(eobj, "npar_psi")

  par_sigma <- par[seq_len(npar_sigma)]
  par_psi  <- par[seq_len(npar_psi) + npar_sigma]

  sigma <- diag(ndim_j)
  sigma[lower.tri(sigma)] <- par_sigma
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]

  psi <- diag(z2r(par_psi))

  sigma0 <-
    matrix(solve(diag(ndim_j^2) - kronecker(psi, psi),
                 c(sigma)),
           ncol = ndim_j)

  ar_blocks <- lapply(0:(ndim_t - 1), function(t) {
    s_ar <- if (t == 0) diag(ndim_j) else psi ^ t
    tcrossprod(sigma0, s_ar)
  })
  S <- toeplitz.block(ar_blocks)
  list(S)
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

  npar1 <- attr(eobj, "npar_sigma")/nlev
  npar1_psi <- attr(eobj, "npar_psi")/nlev
  par_sigma <- par[seq_len(npar_sigma)]
  par_psi <- par[seq_len(npar_psi)+ npar_sigma]

  sigma <- lapply(seq_len(nlev), function(l) {
    S1 <- diag(ndim_j)
    S1[lower.tri(S1)] <- par_sigma[(l - 1) * npar1 + seq_len(npar1)]
    S1[upper.tri(S1)] <- t(S1)[upper.tri(S1)]
    S1
  })

  ar_blocks <- lapply(seq_len(nlev), function(l) {lapply(0:(ndim_t - 1), function(t) {
    if (t == 0) sigma[[l]] else matrix(0, ndim_j, ndim_j)
  })})
  S <- lapply(seq_len(nlev), function(l) toeplitz.block(ar_blocks[[l]]))
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

  npar1 <- attr(eobj, "npar_sigma")/nlev
  npar1_psi <- attr(eobj, "npar_psi")/nlev
  par_sigma <- par[seq_len(npar_sigma)]
  par_psi   <- par[seq_len(npar_psi)+ npar_sigma]

  psi <- lapply(seq_len(nlev), function(l) {
    diag(z2r(par_psi[(l - 1) * npar1_psi + seq_len(npar1_psi)]))
  })

  sigma0 <- lapply(seq_len(nlev), function(l) {
    matrix(solve(diag(ndim_j^2) - kronecker(psi[[l]], psi[[l]]),
                 c(diag(ndim_j))),
           ncol = ndim_j)
  })

  ar_blocks <- lapply(seq_len(nlev), function(l) {lapply(0:(ndim_t - 1), function(t) {
    s_ar <- if (t == 0) diag(ndim_j) else psi[[l]] ^ t
    sigma0[[l]] %*% s_ar
  })})

  S <- lapply(seq_len(nlev), function(l) toeplitz.block(ar_blocks[[l]]))
  S
}


##############################
#### Methods for cor_MMO3 ####
##############################
# eobj <- rho$error.structure
# y <- rho$y
start_values.cor_MMO3 <- start_values.cor_MMO3_cross <- start_values.cor_MMO3_ar1 <-
  function(eobj, y) {
    ## builds starting values for the correlation structure
    rep(1, attr(eobj, "npar"))
}

# eobj <- error.structure
# data <- data.mvord

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

    npar_sigma <-  ndim_j * (ndim_j - 1)/2 * NCOL(attr(eobj, "covariate"))
    npar_psi <-  ndim_j * NCOL(attr(eobj, "covariate"))
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

    npar_sigma <-  ndim_j * (ndim_j - 1)/2 * NCOL(attr(eobj, "covariate"))
    npar_psi <-  0#ndim_j * NCOL(attr(eobj, "covariate"))
    attr(eobj, "npar_sigma") <- npar_sigma
    attr(eobj, "npar_psi") <- npar_psi
    npar <- npar_sigma + npar_psi
    attr(eobj, "npar") <- npar

    r <- eobj$value
    if (is.null(eobj$fixed)) eobj$fixed  <- FALSE
    attr(eobj, "npar.cor") <- ifelse(eobj$fixed, 0, npar)
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

    npar_sigma <-  0 # ndim_j * (ndim_j - 1)/2 * NCOL(attr(eobj, "covariate"))
    npar_psi <-  ndim_j * NCOL(attr(eobj, "covariate"))
    attr(eobj, "npar_sigma") <- npar_sigma
    attr(eobj, "npar_psi") <- npar_psi
    npar <- npar_sigma + npar_psi
    attr(eobj, "npar") <- npar

    r <- eobj$value

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
    ndim_j <- attr(eobj, "ndim_j")
    ndim_t <- attr(eobj, "ndim_t")
    npar_sigma <- attr(eobj, "npar_sigma")
    npar_psi <- attr(eobj, "npar_psi")
    covar <- attr(eobj, "covariate")
    nlev <- NCOL(covar)
    npar1 <- attr(eobj, "npar_sigma")/nlev
    tpar_sigma <- tpar[seq_len(npar_sigma)]
    tpar_psi <- tpar[seq_len(npar_psi)+ npar_sigma]

    Q <- make_Q_cor_MMO3(tpar, eobj)

    # corr_pars <- Q[[1]][lower.tri(Q[[1]])]
#
    corr_pars <- sapply(Q, function(x) {
      x <- cov2cor(x)
      x[lower.tri(x)]
      })
    # if (npar1 == 1) dim(corr_pars) <- c(1, nlev)


    rVec <- tcrossprod(covar, corr_pars)#tcrossprod(covar, s[lower.tri(s)])
    sdQ  <- sqrt(diag(Q[[1]]))
    return(list(rVec = rVec, sdVec = sdQ))
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
    nlev <- NCOL(covar)
    npar1 <- attr(eobj, "npar_sigma")/nlev
    npar1_psi <- attr(eobj, "npar_psi")/nlev
    tpar_sigma <- tpar[seq_len(npar_sigma)]
    tpar_psi <- tpar[seq_len(npar_psi)+ npar_sigma]
    sigma <- lapply(seq_len(nlev), function(l) {
      nu <- tpar_sigma[(l - 1) * npar1 + seq_len(npar1)]
      angles <- pi * exp(nu)/(1 + exp(nu))
      cosmat <- diag(ndim_j)
      cosmat[lower.tri(cosmat)] <- if(length(angles)>0) cos(angles) else 0
      S1 <- matrix(0, nrow = ndim_j, ncol = ndim_j)
      S1[, 1L] <- 1
      S1[lower.tri(S1, diag = T)][-(1:ndim_j)] <- if(length(angles)>0) sin(angles) else 1
      #S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
      tLmat <- sapply(seq_len(ndim_j),
                      function(j) cosmat[j, ] * cumprod(S1[j, ]))
      sigma <- crossprod(tLmat)
      sigma
    })


    psi <- lapply(seq_len(nlev), function(l) {
      tpar_psi1 <- if (npar1_psi==0) rep(0, ndim_j) else tpar_psi[(l - 1) * npar1_psi + seq_len(npar1_psi)]
      diag(z2r(tpar_psi1))
    })

    ## names
    ynames <- attr(eobj, "ynames")
    ind <- combn(ndim_j,2)
    # row_i <- matrix(rep(seq_len(ndim_j), each = ndim_j), ncol = ndim_j, byrow = TRUE)
    # col_i <- matrix(rep(seq_len(ndim_j), each = ndim_j), ncol = ndim_j)
    # ind <-rbind(row_i[lower.tri(row_i)], col_i[lower.tri(col_i)])
    if (eobj$formula == ~1) {
      ## correlation names
      names.sigma <-
        sapply(seq_len(NCOL(ind)), function(j)
          sprintf("corr %s %s", seq_len(ndim_j)[ind[1,j]],
                  seq_len(ndim_j)[ind[2,j]]))
      names.psi <- paste0("psi ",1:ndim_j)
    } else { ## if factor dependent
      ## correlation names
      names.corr.pair <- apply(ind, 2, function(i)
        paste(i, collapse = " "))
      names.sigma <- paste("corr", rep(colnames(covar), each = NCOL(ind)),
                           rep(names.corr.pair, ncol(covar)))
      names.psi <- paste("psi", rep(colnames(covar), each = NCOL(ind)),
                         rep(1:ndim_j, ncol(covar)))
    }
    attr(eobj, "par") <- unlist(c(if (npar1!=0) lapply(sigma, function(x) x[lower.tri(x)]),
                                  if (npar1_psi!=0) lapply(psi, function(x) diag(x))))#corr_vec
    attr(eobj, "parnames") <-   c(if (npar1==0) NULL else names.sigma, if (npar1_psi==0) NULL else names.psi)
  }
  eobj
}
