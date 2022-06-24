# formula = MMO3(rating, firm_id, year, rater) ~ 0 + LR + LEV + PR + RSIZE + BETA
# #formula = MMO(rating, firm_id, rater) ~ 0 + LR + LEV + PR + RSIZE + BETA
#
# index = rho$index
# y.names = rho$response.name
# x.names = names.rest
# y.levels = rho$response.levels
# response.names = rho$response.names

mvord_data_MMO3 <- function(data, index, y.names, x.names,
                       y.levels, response.names) {
  ## check if response is ordered factor. Set response levels accordingly
  if (is.null(y.levels) & is.ordered(data[,y.names])){
    y.levels <- rep(list(levels(data[,y.names])), length(response.names))
  }
  df <- list()
  data.split.y <- split(data[,c(y.names, index[1])],
                        factor(data[, "MMO3_id"], levels = response.names))
  data.split.x <- split(data[, c(x.names, index[1])],
                        factor(data[, "MMO3_id"], levels = response.names))
  #set colnames (otherwise warning due to identical colnames in reduce)
  for (j in seq_along(data.split.y)) {
    colnames(data.split.y[[j]])[1] <- response.names[j]
    colnames(data.split.x[[j]]) <- c(paste(x.names,j, sep = "."), index[1])
  }


  df$y <- Reduce(function(...) merge(..., by = index[1], all = TRUE),
                 data.split.y)
  subject_id_names <- df$y[,index[1]]
  df$y <- df$y[, - match(index[1], colnames(df$y)), drop = FALSE]
  if (is.null(y.levels)) {
    df$y <- cbind.data.frame(lapply(df$y, ordered))
    df$ylevels <- lapply(seq_len(ncol(df$y)), function(j) levels(df$y[,j]))
  } else {
    df$ylevels <- y.levels
    for (j in seq_along(y.levels)) {
      ## check levels for each response
      #if (!all(levels(df$y[, j]) %in% df$ylevels[[j]]))
      if (!all(unique(df$y[, j]) %in% c(NA,df$ylevels[[j]])))
        warning("levels of response do not all match with response levels", call. = FALSE)
      #      if (!all(y.levels[[j]] %in% unique(df$y[, j])))
      #        warning(sprintf("For response %i, not all response
      #          levels are observed. Model might be non-identifiable if
      #          the thresholds for this response are not restricted.", j),
      #        call.=FALSE)
      df$y[, j] <- ordered(df$y[, j], levels = y.levels[[j]])
    }
  }
  rownames(df$y) <- subject_id_names

  xdatadf <- Reduce(function(...) merge(...,by = index[1], all = TRUE), data.split.x)
  rownames(xdatadf) <- subject_id_names
  xdatadf <- xdatadf[, -match(index[1], colnames(xdatadf)), drop = FALSE]
  df$x <- lapply(1:length(response.names), function(i) {
    tmp <- xdatadf[,(i - 1) * length(x.names) + seq_along(x.names), drop = FALSE]
    names(tmp) <- x.names
    tmp
  })
  names(df$x) <- response.names
  df
}


initialize_MMO3 <- function(rho, formula, data, error.structure, contrasts){
  rho$function.name <- "mvord3"
  rho$formula <- as.formula(paste(formula[[2]][[2]], paste(as.character(formula[-2]), collapse = " ")))
  if (length(formula[[2]]) == 2){
    rho$index <- colnames(data)[1:2]
  } else rho$index <- c(as.character(formula[[2]][[3]]), as.character(formula[[2]][[4]]), as.character(formula[[2]][[5]]))
  rho$response.name <- all.vars(rho$formula[[2]])

  tmp <- data[, rho$index[c(2,3)]]

  data$MMO3_id <- paste0(data[, rho$index[2]], data[, rho$index[3]])
  #rho$response.names_tmp <- lapply(2:3, function(x) unique(data[, rho$index[x]]))#apply(data[, rho$index[c(2,3)]],2, unique)#levels(ordered(data[, rho$index[c(2,3)]]))
  rho$response.names_tmp[[1]] <- unique(data[, rho$index[2]])#TODO: check levels(data[, rho$index[2]])#
  rho$response.names_tmp[[2]] <- unique(data[, rho$index[3]])

  attr(error.structure, "ndim_t") <- length(rho$response.names_tmp[[1]])
  attr(error.structure, "ndim_j") <- length(rho$response.names_tmp[[2]])
  rho$response.names <- do.call("c", lapply(1:length(rho$response.names_tmp[[1]]), function(i) paste0(rho$response.names_tmp[[1]][i], rho$response.names_tmp[[2]])))
  ## checks specific to mvord
  if(length(rho$response.name) > 1) stop("only one response needed", call.=FALSE)
  if (any(is.na(data[, rho$index[1]]))) stop("Missing values are not allowed in the subject index.", call.=FALSE)
  if (any(is.na(data[, rho$index[2]]))) stop("Missing values are not allowed in the measurement index.", call.=FALSE)
  if (any(is.na(data[, rho$index[3]]))) stop("Missing values are not allowed in the measurement index.", call.=FALSE)
  if (!is.null(rho$response.levels) & length(rho$response.levels) != length(rho$response.names))
    stop("Length of response levels must be equal to the number of responses.", call.=FALSE)

  check_args_input1(rho, data) ## used also for mvord2

  if(!is.null(rho$weights.name)){
    if(any(is.na(data[,rho$weights.name]))) {
      data[,rho$weights.name][is.na(data[,rho$weights.name])] <- 0
      warning("Weights with values of NA are set to 0.")
    }
  }
  ## check if more than one response --

  rho$intercept <- ifelse(attr(terms.formula(rho$formula),
                               "intercept") == 1, TRUE, FALSE)
  rho$x.names <- c(all.vars(rho$formula[[3]]),
                   all.vars(formula(error.structure)[[2]]))
  if (any(is.na(data[, rho$x.names]))) stop("Missing values in the covariates are not allowed.")

  names.rest <- colnames(data)
  data.mvord <- mvord_data_MMO3(data, index = rho$index, y.names = rho$response.name, x.names = names.rest,
                           y.levels = rho$response.levels,
                           response.names = rho$response.names)

  rho$levels <- data.mvord$ylevels

  rho$y <- data.mvord$y
  rho$y.names <- colnames(rho$y)

  rho$ndim <- ncol(rho$y)
  rho$n <- nrow(rho$y)
  rho$x <- lapply(seq_len(rho$ndim), function(j) {
    rhs.form <- rho$formula
    rhs.form[[2]] <- NULL
    new.rhs.form <- update(rhs.form, ~  . + 1)

    tmp <- suppressWarnings(model.matrix(new.rhs.form,
                                         model.frame(new.rhs.form,  data.mvord$x[[j]], na.action = function(x) x),
                                         contrasts.arg = contrasts))
    attribute <- attr(tmp, "assign")
    if(rho$intercept == FALSE){
      attribute <- attribute[-1]
      tmp <- tmp[,-1, drop = F]
    }
    attr(tmp, "assign") <- attribute
    tmp
  })

  if (is.null(rho$weights.name)) {
    rho$weights <- rep(1, nrow(rho$y))
  } else {
    tmp <- sapply(data.mvord$x, function(j) as.numeric(j[,rho$weights.name]))
    rho$weights <- apply(tmp,1,function(x) unique(x[!is.na(x)]))
    if(is.list(rho$weights)) stop("Weights need to be constant across multiple measurements", call.=FALSE)
    if(any(rho$weights < 0)) stop("Weights must be non-negative", call.=FALSE)
  }
  ## initialize error structure
  rho$error.structure <- init_fun(error.structure, data.mvord, contrasts)
  ## set offset
  if (is.null(rho$offset)) {
    rho$offset <- lapply(seq_len(rho$ndim), function(j) {
      rhs.form <- rho$formula
      rhs.form[[2]] <- NULL
      mf <- model.frame(rhs.form, data.mvord$x[[j]],
                        na.action = function(x) x)
      #mf[is.na(mf)] <- 0
      if(NCOL(mf) > 0){
        for(k in seq_along(NCOL(mf))){
          if(is.numeric(mf[,k]))
            mf[is.na(mf[,k]),k] <- 0
        }
      }
      model.offset(mf)
    })
  }
  check_args_input2(rho, data)
  rho
}
