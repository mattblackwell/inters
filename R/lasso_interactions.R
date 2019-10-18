##' Post-double selection estimator for interactions
##'
##' \code{post_ds_interaction} applies post-double selection to the
##' estimation of an interaction in a linear model. 
##' 
##' @param data data.frame to find the relevant variables.
##' @param treat string with the name of the treatment variable.
##' @param moderator string with the name of the moderating variable.
##' @param outcome string with the name of the outcome variable. 
##' @param control_vars vector of strings with the names of the
##'   control variables to include.
##' @param panel_vars vector of strings with the names of categorical
##'   variables to include as fixed effects.
##' @param moderator_marg logical indicating if the lower-order term
##'   of the moderator should be included ()
##' @param cluster string with the name of the cluster variable. 
##' @return Returns an object of the class \code{lm} with an
##'   additional \code{clustervcv} object containing the
##'   cluster-robust variance matrix estimate when \code{cluster} is 
##'   provided.
##'
##' @details The \code{post_ds_interaction} implements the post-double
##'   selection estimator of Belloni et al (2014) as applied to
##'   interactions, which was proposed by Blackwell and Olson (2019).
##'   Variables passed to \code{panel_vars} are considered factors
##'   for fixed effects and whose "base effects" are removed by
##'   demeaning all variables by those factors. Interactions between
##'   the moderator and all variables (including the factors generated
##'   by \code{panel_vars}) are generated and passed to the
##'   post-double selection procedure. Base terms for the treatment,
##'   moderator, and control variables are forced to be included in
##'   the final post-double selection OLS. The \code{cluster} argument
##'   adjusts the lasso 
##' 
##' 
##' @references Alexandre Belloni, Victor Chernozhukov, Christian
##'   Hansen, Inference on Treatment Effects after Selection among
##'   High-Dimensional Controls, The Review of Economic Studies,
##'   Volume 81, Issue 2, April 2014, Pages 608-650,
##'   \url{https://doi.org/10.1093/restud/rdt044}
##'
##' Matthew Blackwell and Michael Olson. "Reducing Model Misspectation
##'   and Bias in the Estimation of Interactions." Working Paper,
##'   2019. 
##' @export
##' @importFrom stats as.formula model.matrix

post_ds_interaction <- function(data, treat, moderator, outcome, control_vars,
                                panel_vars = NULL, moderator_marg = TRUE,
                                cluster = NULL) {

  # create y and x matrices
  y_mat <- as.matrix(data[, outcome])
  y_mat <- y_mat - mean(y_mat)
  x_mat <- as.matrix(data[, c(control_vars)])
  t_mat <- cbind(data[, treat], data[, treat] * data[, moderator])
  colnames(t_mat) <- c(treat, paste(c(treat, moderator), collapse = "_"))

  n <- nrow(x_mat)
  
  if (!is.null(cluster)) {
    cl <- data[, cluster]
  } else {
    cl <- NULL
  }
  
  if (!is.null(panel_vars)) {
    fnames <- paste("factor(", panel_vars, ")", sep = "")
    contr.list <- list(contr.sum, contr.sum)
    names(contr.list) <- fnames
    panel_form <- as.formula(paste("~", paste(fnames, collapse = " + ")))
    panel_mat <- model.matrix(panel_form, data = data,
                              contrasts.arg = contr.list)[, -1]
    ## x_mat <- cbind(x_mat, panel_mat)
    x_pan_int <- panel_mat * data[, moderator]
    colnames(x_pan_int) <- paste(colnames(panel_mat), moderator, sep = "_")

  }

  ## create interaction matrix
  
  x_int_mat <- x_mat * data[, moderator]
  if (dim(x_int_mat)[2])
    colnames(x_int_mat) <- paste(colnames(x_mat), moderator, sep = "_")
  


  ## add interaction matrix to x matrix
  if (moderator_marg) {
    x_mat <- cbind(data[, moderator], x_mat)
    colnames(x_mat)[1] <- moderator
  }

  if (!is.null(panel_vars)) {
    forced <- c(rep(TRUE, ncol(x_mat)), rep(FALSE, ncol(x_int_mat)),
                rep(FALSE, ncol(x_pan_int)))
    x_mat <- cbind(x_mat, x_int_mat, x_pan_int)
    hold <- lfe::demeanlist(cbind(y_mat, t_mat, x_mat),
                            fl = lapply(data[, panel_vars], factor))
    y_mat <- hold[,1]
    t_mat <- hold[,2:3]
    x_mat <- hold[,-c(1:3)]    
  } else {
    forced <- c(rep(TRUE, ncol(x_mat)), rep(FALSE, ncol(x_int_mat)))
    x_mat <- cbind(x_mat, x_int_mat)
    x_mat <- scale(x_mat, scale = FALSE)
  }
  

  reg_out <- post_double_selection(x = x_mat, d = t_mat, y = y_mat, cl = cl,
                                   panels = panel_mat, forced = forced)
  return(reg_out)
}

##' @importFrom stats lm contr.sum
post_double_selection <- function(x, d, y, cl, panels, forced) {
  k <- dim(d)[2]
  p <- dim(x)[2]
  n <- dim(x)[1]
  selected <- matrix(NA, nrow = p, ncol = k + 1)
  selected[, 1] <- rlasso_cluster(x = x, y = y, cl = cl,
                                  intercept = FALSE)$index
  for (j in 1:k) {
    selected[, j + 1] <- rlasso_cluster(x = x, y = d[, j], cl = cl,
                                   intercept = FALSE)$index
  }
  selected <- (rowSums(selected) + forced) > 0

  selected_x_df <- as.data.frame(cbind(d, x[, selected]))
  reg_out <- lm(y ~ ., data = selected_x_df)
  if (!is.null(cl)) {
    kd <- sum(selected) + k
    if (!is.null(panels)) {
      K <- kd + ncol(panels)
    } else {
      K <- kd
    }
    cls <- unique(cl)
    ncl <- length(cls)
    bread <- summary(reg_out)$cov.unscaled * n
    meat <- matrix(0, nrow = kd, ncol = kd)
    uj <- apply(sandwich::estfun(reg_out), 2,
                function(x) tapply(x, cl, sum))
    meat <- crossprod(uj) / n
    dfc <- (ncl / (ncl - 1)) * (n - 1) / (n - K)
    vcv <- (1 / n) * dfc * bread %*% meat %*% bread
    rownames(vcv) <- colnames(vcv) <- colnames(bread)
    reg_out$clustervcv <- vcv
  }
  return(reg_out)
}

##' @importFrom stats lm.fit qnorm sd var
rlasso_cluster <- function(x, y, cl, intercept = TRUE, c_val = 1.1, gamma,
                           num_iter = 15, tol = 10^-5) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- dim(x)[1]
  p <- dim(x)[2]

  if (missing(gamma)) gamma <- 0.1 / log(n)

  ymean <- mean(y)
  if (intercept) {
    xmeans <- colMeans(x)
    x <- scale(x, center = xmeans, scale = FALSE)
    y <- y - ymean
  }
  eps0 <- lm.fit(x, y)$residuals
  lambda0 <- 2 * c_val * sqrt(n) * qnorm(1 - gamma / (2 * p * 1))
  if (is.null(cl)) {
    phi0 <- phi1 <- sqrt(1/n * t(t(eps0^2) %*% x^2))
  } else {
    uj <- apply(eps0 * x, 2, function(x) tapply(x, cl, sum))
    phi0 <- phi1 <- sqrt(diag(crossprod(uj) / n))
  }

  pen_norm <- sum(phi0) / p
  count <- 1
  s0 <- sd(y)
  s1 <- Inf
  while (count <= num_iter & (abs(s0 - s1) > tol)) {
    if (count == 1) {
      cf_tmp <- glmnet::glmnet(x = x, y = y, lambda = pen_norm * lambda0/(4 * n),
                               penalty.factor = phi1, intercept = FALSE,
                               standardize = FALSE)$beta
    } else {
      s0 <- s1
      cf_tmp <- glmnet::glmnet(x = x, y = y, lambda = pen_norm * lambda0/(2 * n),
                               penalty.factor = phi1, intercept = FALSE,
                               standardize = FALSE)$beta
    }
    cf_tmp <- as.numeric(cf_tmp)
    cf_tmp[is.na(cf_tmp)] <- 0
    s_vars <- (abs(cf_tmp) > 0)
    x1 <- x[, s_vars, drop = FALSE]

    if (ncol(x1) == 0) {
      if (intercept) {
        cf <- c(ymean, rep(0, p))
        names(cf) <- c("(Intercept)", colnames(x))
      } else {
        cf <- rep(0, p)
        names(cf) <- colnames(x)
      }
      out <- list(coefficients = cf, index = rep(FALSE, p), lambda = phi1 * lambda0,
                  lambda0 = lambda0,
                  loadings = phi1, iter = count, residuals = y - ymean, sigma = var(y),
                  c_val = c_val, gamma = gamma)
      out$tss <- out$rss <- sum((y - ymean)^2)
      out$dev <- y - ymean
      return(out)
    }

    ehat <- lm.fit(x = x1, y = y)$residuals
    s1 <- sd(ehat)

    ## update pentalties
    if (is.null(cl)) {
      Ups1 <- sqrt(1/n * t(t(ehat^2) %*% x^2))
    } else {
      uj <- apply(ehat * x, 2, function(x) tapply(x, cl, sum))
      Ups1 <- sqrt(diag(crossprod(uj) / n))
    }
    count <- count + 1
  }

  names(cf_tmp) <- names(s_vars) <- colnames(x)
  if (intercept) {
    int_value <- ymean - sum(xmeans * cf_tmp)
    cf <- c(int_value, cf_tmp)
    names(cf)[1] <- "(Intercept)"
  } else {
    cf <- cf_tmp
  }
  
  out <- list(coefficients = cf, index = s_vars, lambda = Ups1 * lambda0, lambda0 = lambda0,
              loadings = Ups1, iter = count, residuals = as.vector(ehat), sigma = s1,
              c_val = c_val, gamma = gamma)
  out$tss <- sum((y - ymean)^2)
  out$rss <- sum(ehat^2)
  out$dev <- y - ymean
  return(out)
}

