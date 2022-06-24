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
#' @param blocks a list of matrices
#'              Sigma, Psi * Sigma, Psi^2 * Sigma, ...,Psi^(T-1) * Sigma
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
  nlev <- NCOL(covar)
  npar1 <- attr(eobj, "npar_sigma")/nlev
  npar1_psi <- attr(eobj, "npar_psi")/nlev
  tpar_sigma <- tpar[seq_len(npar_sigma)]
  tpar_psi <- tpar[seq_len(npar_psi)+ npar_sigma]

  sigma <- lapply(seq_len(nlev), function(l) {
    nu <- tpar_sigma[(l - 1) * npar1 + seq_len(npar1)]
    angles <- pi * exp(nu)/(1 + exp(nu))
    cosmat <- diag(ndim_j)
    cosmat[lower.tri(cosmat)] <- cos(angles)
    S1 <- matrix(0, nrow = ndim_j, ncol = ndim_j)
    S1[, 1L] <- 1
    S1[lower.tri(S1, diag = T)][-(1:ndim_j)] <- sin(angles)
    #S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
    tLmat <- sapply(seq_len(ndim_j),
                    function(j) cosmat[j, ] * cumprod(S1[j, ]))
    sigma <- crossprod(tLmat)
    sigma
  })
  psi <- lapply(seq_len(nlev), function(l) {diag(z2r(tpar_psi[(l - 1) * npar1_psi + seq_len(npar1_psi)]))})

  ar_blocks <- lapply(seq_len(nlev), function(l) {lapply(0:(ndim_t - 1), function(t) {
    s_ar <- if (t == 0) diag(ndim_j) else psi[[l]] ^ t
    crossprod(s_ar, sigma[[l]])
  })})
  S <- lapply(seq_len(nlev), function(l) toeplitz.block(ar_blocks[[l]]))
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

  psi <- lapply(seq_len(nlev), function(l) {diag(z2r(par_psi[(l - 1) * npar1_psi + seq_len(npar1_psi)]))})

  ar_blocks <- lapply(seq_len(nlev), function(l) {lapply(0:(ndim_t - 1), function(t) {
    s_ar <- if (t == 0) diag(ndim_j) else psi[[l]] ^ t
    crossprod(s_ar, sigma[[l]])
  })})
  S <- lapply(seq_len(nlev), function(l) toeplitz.block(ar_blocks[[l]]))
  S
}


#################################
#### Methods for cor_general ####
#################################
# eobj <- rho$error.structure
# y <- rho$y
start_values.cor_MMO3 <- function(eobj, y) {
  ## builds starting values for the correlation structure
  tmp <- rep(1, attr(eobj, "npar"))

  ndim_j <- attr(eobj, "ndim_j")
  covar <- attr(eobj, "covariate")
  nlev <- NCOL(covar)

  # ### starting values with polychoric correlations
  # sigma_tmp <- diag(ndim_j)
  # for(i in 1:(ndim_j-1)){
  #   for(j in (i+1):ndim_j){
  #     sigma_tmp[i,j] <- suppressWarnings(polycor::polychor(y[,i], y[,j]))
  #   }
  # }
  # sigma_tmp[lower.tri(sigma_tmp)] <- t(sigma_tmp)[lower.tri(sigma_tmp)]
  # sigma_tmp[is.na(sigma_tmp)] <- 0.1
  # if(!all(eigen(sigma_tmp)$values >0)){
  #   sigma_tmp <- suppressWarnings(Matrix::nearPD(sigma_tmp)$mat)
  # }
  # tsigma <- suppressWarnings(backtransf_sigmas(sigma_tmp))
  # tsigma[is.na(tsigma)] <- 0
  # psi_tmp <- sapply(seq_len(ndim_j), function(j) suppressWarnings(polycor::polychor(y[,j], y[,j + ndim_j])))
  # tpsi <- atanh(psi_tmp)
  # tmp <- c(rep(tsigma, nlev), rep(tpsi, nlev))
  # ## TODO for the given values
  tmp
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

build_error_struct.cor_MMO3 <-
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


    corr_pars <- sapply(Q, function(x) x[lower.tri(x)])

   # if (npar1 == 1) dim(corr_pars) <- c(1, nlev)


    rVec <- tcrossprod(covar, corr_pars)#tcrossprod(covar, s[lower.tri(s)])#tcrossprod(covar, corr_pars)
    sd <- rep(1, ndim)
    return(list(rVec = rVec, sdVec = sd))
  }

finalize_fun.cor_MMO3 <-
  function(eobj, tpar)
  {
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
        cosmat[lower.tri(cosmat)] <- cos(angles)
        S1 <- matrix(0, nrow = ndim_j, ncol = ndim_j)
        S1[, 1L] <- 1
        S1[lower.tri(S1, diag = T)][-(1:ndim_j)] <- sin(angles)
        #S1[-1L, -1L][lower.tri(S1[-1L, -1L], diag = T)] <- sin(angles)
        tLmat <- sapply(seq_len(ndim_j),
                        function(j) cosmat[j, ] * cumprod(S1[j, ]))
        sigma <- crossprod(tLmat)
        sigma
      })


      psi <- lapply(seq_len(nlev), function(l) {diag(z2r(tpar_psi[(l - 1) * npar1_psi + seq_len(npar1_psi)]))})

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
      attr(eobj, "par") <- unlist(c(lapply(sigma, function(x) x[lower.tri(x)]), lapply(psi, function(x) diag(x))))#corr_vec
      attr(eobj, "parnames") <-   c(names.sigma, names.psi)
    }
    eobj
  }
