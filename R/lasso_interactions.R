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
##' @param method string indicating which method to use. The default
##' is \code{"double selection"} selects variables based on the
##' outcome and treatment/interaction variables and \code{"single
##' selection"} only selects on the outcome.
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
##' @examples
##' data(remit)
##'
##' controls <- c("l1gdp", "l1pop", "l1nbr5", "l12gr", "l1migr",
##' "elec3")
##' 
##' post_ds_out <- post_ds_interaction(
##'   data = remit, treat = "remit",
##'   moderator = "dict", outcome = "Protest",
##'   control_vars = controls,
##'   cluster = "caseid"
##' )
##'
##' @references Alexandre Belloni, Victor Chernozhukov, Christian
##'   Hansen, Inference on Treatment Effects after Selection among
##'   High-Dimensional Controls, The Review of Economic Studies,
##'   Volume 81, Issue 2, April 2014, Pages 608-650,
##'   \doi{10.1093/restud/rdt044}
##'
##' Matthew Blackwell and Michael Olson.. "Reducing Model Misspectation
##'   and Bias in the Estimation of Interactions." Political Analysis,
##' 2021. 
##' @export
##' @importFrom stats as.formula model.matrix coef cor hatvalues
##' residuals model.frame reformulate

post_ds_interaction <- function(data, treat, moderator, outcome, control_vars,
                                panel_vars = NULL, moderator_marg = TRUE,
                                cluster = NULL, method = "double selection") {


  if (!is.null(panel_vars)) {
    if (any(sapply(data[panel_vars], is.numeric))) {
      num_pan <- which(sapply(data[panel_vars], is.numeric))
      panel_vars[num_pan] <- paste("factor(", panel_vars[num_pan], ")", sep = "")
    }
  } 

  all_vars <- c(treat, moderator, control_vars, panel_vars, cluster)
  formula <- stats::reformulate(all_vars, response = outcome)
  mf <- stats::model.frame(formula, data)
  # create y and x matrices
  
  y_mat <- as.matrix(mf[[outcome]])
  y_mat <- y_mat - mean(y_mat)
  control_form <- stats::reformulate(control_vars, intercept = FALSE)
  x_mat <- model.matrix(control_form, data = mf)
  t_mat <- cbind(mf[[treat]], mf[[treat]] * mf[[moderator]])
  colnames(t_mat) <- c(treat, paste(c(treat, moderator), collapse = "_"))


  if (!is.null(cluster)) {
    cl <- mf[[cluster]]
  } else {
    cl <- NULL
  }

  if (!is.null(panel_vars)) {
    contr.list <- list(contr.sum, contr.sum)
    names(contr.list) <- panel_vars
    panel_form <- stats::reformulate(panel_vars, intercept = FALSE)
    panel_mat <- model.matrix(panel_form, data = mf,
                              contrasts.arg = contr.list)
  ## x_mat <- cbind(x_mat, panel_mat)
    x_pan_int <- panel_mat * mf[[moderator]]
    colnames(x_pan_int) <- paste(colnames(panel_mat), moderator, sep = "_")

  } else {
    panel_mat <- NULL
  }

  ## create interaction matrix
  x_int_mat <- x_mat * mf[[moderator]]
  if (dim(x_int_mat)[2])
    colnames(x_int_mat) <- paste(colnames(x_mat), moderator, sep = "_")


  ## add interaction matrix to x matrix
  if (moderator_marg | method == "partialing out") {
    x_mat <- cbind(mf[[moderator]], x_mat)
    colnames(x_mat)[1] <- moderator
  }

  if (!is.null(panel_vars)) {
    forced <- c(rep(TRUE, ncol(x_mat)), rep(FALSE, ncol(x_int_mat)),
                rep(FALSE, ncol(x_pan_int)))
    x_mat <- cbind(x_mat, x_int_mat, x_pan_int)
    hold <- fixest::demean(cbind(y_mat, t_mat, x_mat),
                           f = mf[panel_vars])
    y_mat <- hold[, 1]
    t_mat <- hold[, 2:3]
    x_mat <- hold[, -c(1:3)]
  } else {
    forced <- c(rep(TRUE, ncol(x_mat)), rep(FALSE, ncol(x_int_mat)))
    x_mat <- cbind(x_mat, x_int_mat)
    x_mat <- scale(x_mat, scale = FALSE)
  }

  reg_out <- post_double_selection(x = x_mat, d = t_mat, y = y_mat,
                                   v = mf[[moderator]], cl = cl,
                                   panels = panel_mat, forced = forced,
                                   method = method)
  return(reg_out)
}

##' @importFrom stats lm contr.sum
post_double_selection <- function(x, d, y, v, cl, panels, forced, method) {
  k <- dim(d)[2]
  p <- dim(x)[2]
  n <- dim(x)[1]
  dd <- scale(d, scale = FALSE)
  if (method %in% c("double selection", "single selection")) {
    if (method == "double selection") {
      selected <- matrix(NA, nrow = p, ncol = k + 1)
      selected[, 1] <- rlasso_cluster(x = cbind(dd, x), y = y, cl = cl,
                                      intercept = FALSE)$index[(k + 1):(k + p)]

      for (j in 1:k) {
        selected[, j + 1] <- rlasso_cluster(x = cbind(dd[, -j], x),
                                            y = dd[, j], cl = cl,
                                            intercept = FALSE)$index[k:(k + p - 1)]
      }
      selected <- cbind(selected, forced)
      colnames(selected) <- c("y", colnames(d), "forced")

    } else {
      selected <- matrix(NA, nrow = p, ncol = 1)
      selected[, 1] <- rlasso_cluster(x = cbind(dd, x), y = y, cl = cl,
                                      intercept = FALSE)$index[(k + 1):(k + p)]
      selected <- cbind(selected, forced)
      colnames(selected) <- c("y", "forced")

    }

    rownames(selected) <- colnames(x)
    selected_all <- (rowSums(selected) > 0)

    selected_x_df <- as.data.frame(cbind(d, x[, selected_all]))
    reg_out <- lm(y ~ ., data = selected_x_df)
    reg_out$selection_matrix <- selected
    mm <- model.matrix(reg_out)
    dd <- mm[, 2:(k + 1)]
    mm <- mm[, setdiff(seq_len(ncol(mm)), 2:(k + 1))]
    res_d <- residuals(lm(dd ~ mm - 1))

    res_y <- residuals(lm(y ~ mm - 1))
    res_y <- res_y - res_d %*% coef(reg_out)[2:(k + 1)]
    p_mat <- res_d %*% solve(crossprod(res_d)) %*% t(res_d)
    res_y <- c(res_y) / (1 - diag(p_mat))
    ef <- res_d * res_y
    bread <- solve(crossprod(res_d) / n)
    tot_p <- k + sum(selected)
  } else if (method == "partialing out") {
    ## with panel models we force V to be here, so we need to drop it
    ## if doesn't vary after demeaning
    ## first get rid of the d/v
    s1 <- sample(1:n, size = round(n / 2))
    s2 <- setdiff(1:n, s1)

    reg_y_s1 <- rlasso_cluster(x = x[s1, ], y = y[s1], cl = cl[s1], intercept = TRUE)
    reg_y_s2 <- rlasso_cluster(x = x[s2, ], y = y[s2], cl = cl[s2], intercept = TRUE)
    r_y_s1 <- y[s1] - cbind(1, x[s1, ]) %*% reg_y_s2$coefficients
    r_y_s2 <- y[s2] - cbind(1, x[s2, ]) %*% reg_y_s1$coefficients

    reg_d_s1 <- rlasso_cluster(x = x[s1, ], y = d[s1, 1], cl = cl[s1], intercept = TRUE)
    reg_d_s2 <- rlasso_cluster(x = x[s2, ], y = d[s2, 1], cl = cl[s2], intercept = TRUE)
    r_d_s1 <- d[s1, 1] - cbind(1, x[s1, ]) %*% reg_d_s2$coefficients
    r_d_s2 <- d[s2, 1] - cbind(1, x[s2, ]) %*% reg_d_s1$coefficients

    reg_dv_s1 <- rlasso_cluster(x = x[s1, ], y = d[s1, 2], cl = cl[s1], intercept = TRUE)
    reg_dv_s2 <- rlasso_cluster(x = x[s2, ], y = d[s2, 2], cl = cl[s2], intercept = TRUE)
    r_dv_s1 <- d[s1, 2] - cbind(1, x[s1, ]) %*% reg_dv_s2$coefficients
    r_dv_s2 <- d[s2, 2] - cbind(1, x[s2, ]) %*% reg_dv_s1$coefficients

    reg_data1 <- data.frame(y = r_y_s1, d = r_d_s1, d_V = r_dv_s1)
    reg_data2 <- data.frame(y = r_y_s2, d = r_d_s2, d_V = r_dv_s2)
    reg_out <- lm(y ~ . - 1, data = reg_data1)
    reg_out2 <- lm(y ~ . - 1, data = reg_data2)

    reg_out$coefficients <- (reg_out$coefficients + reg_out2$coefficients) / 2
    num_sel <- sum((reg_y_s1$index + reg_d_s1$index + reg_dv_s1$index) > 0) + 1
    mm <- model.matrix(reg_out)
    res_y <- residuals(reg_out) / (1 - hatvalues(reg_out))
    ef <- mm * res_y
    bread <- solve(crossprod(mm) / n)
    tot_p <- k + num_sel
  }

  if (!is.null(panels)) tot_p <- tot_p + ncol(panels)
  if (!is.null(cl)) {
    cls <- unique(cl)
    ncl <- length(cls)
    uj <- apply(ef, 2, function(x) tapply(x, cl, sum))
    meat <- crossprod(uj) / n
    dfc <- (ncl / (ncl - 1)) * (n - 1) / (n - tot_p)
  } else {
    meat <- crossprod(ef) / n
    dfc <-  (n / (n - tot_p)) * ((n - 1) / n)
  }
  vcv <- (1 / n) * dfc * bread %*% meat %*% bread
  rownames(vcv) <- colnames(vcv) <- colnames(bread)
  reg_out$vcv <- vcv
  reg_out$se <- sqrt(diag(vcv))
  return(reg_out)
}

##' @importFrom stats lm.fit qnorm sd var
rlasso_cluster <- function(x, y, cl, intercept = TRUE, c_val = 1.1, gamma,
                           num_iter = 15, tol = 10^-5, post = TRUE) {
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
  corr <- apply(x, 2, function(xx) ifelse(sd(xx) > 0, cor(y, xx), 0))
  index <- order(corr, decreasing = TRUE)[1:min(5, p)]
  eps0 <- lm.fit(x[, index, drop = FALSE], y)$residuals
  lambda0 <- 2 * c_val * sqrt(n) * qnorm(1 - gamma / (2 * p * 1))
  if (is.null(cl)) {
    phi0 <- phi1 <- sqrt(1/n * t(t(eps0 ^ 2) %*% x ^ 2))
  } else {
    uj <- apply(eps0 * x, 2, function(x) tapply(x, cl, sum))
    phi0 <- sqrt(diag(crossprod(uj) / n))
  }

  pen_norm <- sum(phi0) / p
  count <- 1
  phi1 <- rep(Inf, p)
  while (count <= num_iter & (max(abs(phi0 - phi1)) > tol)) {
    if (count == 1) {
      cf_tmp <- glmnet::glmnet(x = x, y = y, lambda = pen_norm * lambda0/(4 * n),
                               penalty.factor = phi0, intercept = FALSE,
                               standardize = FALSE)$beta
    } else {
      phi0 <- phi1
      pen_norm <- sum(phi1) / p
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
      phi1 <- sqrt(1 / n * t(t(ehat ^ 2) %*% x ^ 2))
    } else {
      uj <- apply(ehat * x, 2, function(x) tapply(x, cl, sum))
      phi1 <- sqrt(diag(crossprod(uj) / n))
    }
    count <- count + 1
  }

  names(cf_tmp) <- names(s_vars) <- colnames(x)
  if (post) {
    reg <- lm.fit(x = x1, y = y)
    ehat <- reg$residuals
    cf_p <- reg$coefficients
    cf_p[is.na(cf_p)] <- 0
    cf_tmp[s_vars] <- cf_p
  } else {
    ehat <- y - x1 %*% cf_tmp[s_vars]
  }
  cf <- cf_tmp
  s1 <- sd(ehat)

  if (intercept) {
    int_value <- ymean - sum(xmeans * cf)
    cf <- c(int_value, cf)
    names(cf)[1] <- "(Intercept)"
  }

  out <- list(coefficients = cf, index = s_vars, lambda = phi1 * lambda0, lambda0 = lambda0,
              loadings = phi1, iter = count, residuals = as.vector(ehat), sigma = s1,
              c_val = c_val, gamma = gamma)
  out$tss <- sum((y - ymean)^2)
  out$rss <- sum(ehat^2)
  out$dev <- y - ymean
  return(out)
}
