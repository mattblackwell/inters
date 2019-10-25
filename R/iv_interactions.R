##' Estimates principal stratum-specific effects and interactions in a
##' 2x2 factorial experiment
##'
##' This function estimates treatment effects for 2x2 factorial
##' experiments in the face of noncompliance on both factors. A
##' monotonicity assumption is assumed for both treatment-instrument
##' pairs, along with treatment exclusion. See Blackwell (2017) for
##' more details on those assumptions.
##'
##' The procedure uses iterative generalized method of moments (GMM)
##' to  estimate both the proportions of each compliance class (also
##' known as principal strata) and the average potential outcomes
##' within those classes. It also provides estimates of several
##' one-way, joint, and interactive treatment effects within these
##' classes. 
##'
##' Under the above assumptions, the compliance classes are the
##' product of the compliance classes for each treatment-instrument
##' pair. For instance, \code{"cc"} is the class that would comply
##' with both treatments, \code{"ca"} is the class that would comply
##' with the first treatment and always take the second treatment, and
##' \code{"cn"} is the class that would comply with the first
##' treatment and never take the second treatment. In the estimated
##' effects, some effects are defined for the \code{"c-"} or
##' \code{"-c"} compliance classes, which are effects for compliers on
##' one treatment, averaging over the compliance status of the other
##' treatment. Finally, note that treatment effects are only
##' well-defined for compliance classes for which there is compliance
##' on at least one treatment. 
##' 
##' @title IV Estimation of 2x2 Factorial Design
##' @param y vector of outcomes.
##' @param d 2-column matrix of binary treatments.
##' @param z 2-column matrix of binary instruments.
##' @param max_iter maximum number of iterations to conduct when
##'   estimating the GMM. 
##' @param tol criteria for stopping the iterative GMM procedure. 
##' @return A list of class \code{iv2x2} that contains the following
##'   components: 
##' \item{rho}{vector of estimated compliance class
##'   probabilities.}
##' \item{psi}{vector of the estimated conditional mean of the outcome
##'   within the compliance classes.}
##' \item{init}{list of two vectors, with the initial starting values
##'   for the GMM estimation based on a just-identified estimation.}
##' \item{vcov}{estimated asymptotic variance matrix of the combined
##'   \code{rho} and \code{psi} parameters.}
##' \item{tau}{vector of estimated causal effects within levels of the
##'   compliance classes.}
##' \item{tau_se}{vector of estimated standard errors for the
##'   estimated effects in \code{tau}.}
##' @author Matt Blackwell
##' @references Matthew Blackwell (2017) Instrumental Variable Methods
##'   for Conditional Effects and Causal Interaction in Voter
##'   Mobilization Experiments, Journal of the American Statistical
##'   Association, 112:518, 590-599,
##'   \url{https://doi.org/10.1080/01621459.2016.1246363}
##' @export
##' @importFrom stats optim cor 

iv_2x2 <- function(y, d, z, max_iter = 15, tol = 10e-7) {
  if (dim(d)[2] != 2) stop("d does not have two columns")
  if (dim(z)[2] != 2) stop("z does not have two columns")
  
  pstr <- expand.grid(c("c", "a", "n"), c("c", "a", "n"))
  pstr_names <- paste0(pstr$Var1, pstr$Var2)[-9]
  N <- length(y)
  inits <- iv_init(y, d, z)
  dat <- cbind(d, z, y)
  g_inits <- iv_int_moments(theta = unlist(inits), x = dat)
  W0 <- diag(nrow = ncol(g_inits))
  first <- optim(par = unlist(inits), fn = iv_int_loss, gr = iv_int_Qgrad,
                 x = dat, W = W0, method = "BFGS")
  res <- first$par
  delt <- 1000
  count <- 0
  while (count <= max_iter & delt > tol) {
    ghat <- iv_int_moments(res, x = dat)
    Lhat <- (1/N) * t(ghat) %*% ghat
    second <- optim(par = unlist(inits), fn = iv_int_loss, gr = iv_int_Qgrad,
                    x = dat, W = solve(Lhat))
    delt <- sum(abs(res - second$par) / res)
    res <- second$par
    count <- count + 1
  }
  if (count == max_iter & delt > tol) {
    warning("Reached maximum iterations without converging...\n")
  }
  Ghat <- (1/N) * iv_int_grad(res, x = dat)
  bads <- which(abs(cor(Ghat)) > 0.99 & abs(cor(Ghat)) < 1, arr.ind = TRUE)
  if (length(bads)) {
    bad_rhos <- bads[bads[,1] %in% 1:8, 1]
    drops <- grepl(substr(names(res)[bad_rhos], 5, 7), names(res))
    Ghat <- Ghat[, !drops]
  } else {
    drops <- rep_len(FALSE, length(res))
  }
  res_var <- solve(t(Ghat) %*% solve(Lhat) %*% Ghat) / N
  out <- list()
  out$rho <- res[1:8]
  names(out$rho) <- pstr_names
  out$psi <- res[9:24]
  names(out$psi) <- c("11.cc", "11.ac", "11.ca", "11.aa",
                      "10.cc", "10.ac", "10.cn", "10.an",
                      "01.cc", "01.nc", "01.na", "01.ca",
                      "00.cc", "00.nc", "00.cn", "00.nn")
  colnames(res_var) <- c(names(out$rho), names(out$psi))[!drops]
  rownames(res_var) <- c(names(out$rho), names(out$psi))[!drops]
  out$vcov <- res_var
  out$init <- inits
  
  effs <- psi_to_tau(out$psi, out$rho)
  out$tau <- effs$tau
  out$tau_se <- sqrt(diag(tcrossprod(effs$Chat[,!drops] %*% out$vcov, effs$Chat[,!drops])))
  out$tau_se[out$tau_se == 0] <- NA
  class(out) <- "iv2x2"
  return(out)
}


iv_int_moments <- function(theta, x, W) {
  d1 <- x[, 1]
  d2 <- x[, 2]
  z1 <- x[, 3]
  z2 <- x[, 4]
  y <- x[, 5]
  N <- length(y)
  z11 <- z1 * z2
  z10 <- z1 * (1 - z2)
  z01 <- (1 - z1) * z2
  z00 <- (1 - z1) * (1 - z2)
  d11 <- d1 * d2
  d10 <- d1 * (1 - d2)
  d01 <- (1 - d1) * d2
  d00 <- (1 - d1) * (1 - d2)

  pstr <- expand.grid(c("c", "a", "n"), c("c", "a", "n"))
  pstr_names <- paste0(pstr$Var1, pstr$Var2)
  rho <- c(theta[1:8], 1 - sum(theta[1:8]))
  names(rho) <- pstr_names
  psi <- theta[9:24]
  names(psi) <- c("11.cc", "11.ac", "11.ca", "11.aa",
                  "10.cc", "10.ac", "10.cn", "10.an",
                  "01.cc", "01.nc", "01.na", "01.ca",
                  "00.cc", "00.nc", "00.cn", "00.nn")

  ## treatment momemnt conditions
  fvecs <- matrix(NA, nrow = N, ncol = 12)
  colnames(fvecs) <- c("f11.11", "f01.11", "f10.11",
                       "f11.01", "f01.01", "f10.01",
                       "f11.10", "f01.10", "f10.10",
                       "f11.00", "f01.00", "f10.00")

  fvecs[, "f11.11"] <- z11 * (d11 - rho["cc"] - rho["aa"] - rho["ac"] -
                                rho["ca"])
  fvecs[, "f10.11"] <- z11 * (d10 - rho["cn"] - rho["an"])
  fvecs[, "f01.11"] <- z11 * (d01 - rho["nc"] - rho["na"])
  fvecs[, "f11.01"] <- z01 * (d11 - rho["aa"] - rho["ac"])
  fvecs[, "f10.01"] <- z01 * (d10 - rho["an"])
  fvecs[, "f01.01"] <- z01 * (d01 - rho["cc"] - rho["nc"] - rho["na"] -
                                rho["ca"])
  fvecs[, "f11.10"] <- z10 * (d11 - rho["aa"] - rho["ca"])
  fvecs[, "f10.10"] <- z10 * (d10 - rho["cc"] - rho["ac"] - rho["cn"] -
                                rho["an"])
  fvecs[, "f01.10"] <- z10 * (d01 - rho["na"])
  fvecs[, "f11.00"] <- z00 * (d11 - rho["aa"])
  fvecs[, "f10.00"] <- z00 * (d10 - rho["ac"] - rho["an"])
  fvecs[, "f01.00"] <- z00 * (d01 - rho["ca"] - rho["na"])

  ## outcome moment equations
  svecs <- matrix(NA, nrow = N, ncol = 16)
  colnames(svecs) <- c("s11.11", "s01.11", "s10.11", "s00.11",
                       "s11.01", "s01.01", "s10.01", "s00.01",
                       "s11.10", "s01.10", "s10.10", "s00.10",
                       "s11.00", "s01.00", "s10.00", "s00.00")

  svecs[, "s11.11"] <- z11 * (y * d11 - psi["11.cc"] * rho["cc"] -
                                psi["11.ca"] * rho["ca"] -
                                psi["11.ac"] * rho["ac"] -
                                psi["11.aa"] * rho["aa"])
  svecs[, "s10.11"] <- z11 * (y * d10 - psi["10.cn"] * rho["cn"] -
                                psi["10.an"] * rho["an"])
  svecs[, "s01.11"] <- z11 * (y * d01 - psi["01.nc"] * rho["nc"] -
                                psi["01.na"] * rho["na"])
  svecs[, "s00.11"] <- z11 * (y * d00 - psi["00.nn"] * rho["nn"])
  svecs[, "s11.10"] <- z10 * (y * d11 - psi["11.aa"] * rho["aa"] -
                                psi["11.ca"] * rho["ca"])
  svecs[, "s10.10"] <- z10 * (y * d10 - psi["10.cc"] * rho["cc"] -
                                psi["10.ac"] * rho["ac"] -
                                psi["10.cn"] * rho["cn"] -
                                psi["10.an"] * rho["an"])
  svecs[, "s01.10"] <- z10 * (y * d01 - psi["01.na"] * rho["na"])
  svecs[, "s00.10"] <- z10 * (y * d00 - psi["00.nn"] * rho["nn"] -
                                psi["00.nc"] * rho["nc"])
  svecs[, "s11.01"] <- z01 * (y * d11 - psi["11.aa"] * rho["aa"] -
                                psi["11.ac"] * rho["ac"])
  svecs[, "s10.01"] <- z01 * (y * d10 - psi["10.an"] * rho["an"])
  svecs[, "s01.01"] <- z01 * (y * d01 - psi["01.cc"] * rho["cc"] -
                                psi["01.nc"] * rho["nc"] -
                                psi["01.na"] * rho["na"] -
                                psi["01.ca"] * rho["ca"])
  svecs[, "s00.01"] <- z01 * (y * d00 - psi["00.cn"] * rho["cn"] -
                                psi["00.nn"] * rho["nn"])

  svecs[, "s11.00"] <- z00 * (y * d11 - psi["11.aa"] * rho["aa"])
  svecs[, "s10.00"] <- z00 * (y * d10 - psi["10.an"] * rho["an"] -
                                psi["10.ac"] * rho["ac"])
  svecs[, "s01.00"] <- z00 * (y * d01 - psi["01.na"] * rho["na"] -
                                psi["01.ca"] * rho["ca"])
  svecs[, "s00.00"] <- z00 * (y * d00 - psi["00.cc"] * rho["cc"] -
                                psi["00.nn"] * rho["nn"] -
                                psi["00.cn"] * rho["cn"] -
                                psi["00.nc"] * rho["nc"])

  out <- cbind(fvecs, svecs)
  out
}

iv_int_loss <- function(theta, x, W) {
  gmats <- iv_int_moments(theta, x)
  ghats <- colMeans(gmats)
  loss <- crossprod(ghats, W %*% ghats)
  return(loss)
}



psi_to_tau <- function(psi, rho) {

  tau <- rep_len(NA, 14)
  names(tau) <- c("10.00.cc", "10.00.cn", "10.00.c-",
                  "01.00.cc", "01.00.nc", "01.00.-c",
                  "11.01.cc", "11.01.ca", "11.01.c-",
                  "11.10.cc", "11.10.ac", "11.10.-c",
                  "11.00.cc", "int.cc")
  tau["10.00.cc"] <- psi["10.cc"] - psi["00.cc"]
  tau["10.00.cn"] <- psi["10.cn"] - psi["00.cn"]
  tau["10.00.c-"] <- (rho["cc"]/(rho["cc"] + rho["cn"])) * tau["10.00.cc"] +
    (rho["cn"]/(rho["cc"] + rho["cn"])) * tau["10.00.cn"]
  tau["01.00.cc"] <- psi["01.cc"] - psi["00.cc"]
  tau["01.00.nc"] <- psi["01.nc"] - psi["00.nc"]
  tau["01.00.-c"] <- (rho["cc"]/(rho["cc"] + rho["nc"])) * tau["01.00.cc"] +
    (rho["nc"]/(rho["cc"] + rho["nc"])) * tau["01.00.nc"]
  tau["11.01.cc"] <- psi["11.cc"] - psi["01.cc"]
  tau["11.01.ca"] <- psi["11.ca"] - psi["01.ca"]
  tau["11.01.c-"] <- (rho["cc"]/(rho["cc"] + rho["ca"])) * tau["11.01.cc"] +
    (rho["ca"]/(rho["cc"] + rho["ca"])) * tau["11.01.ca"]
  tau["11.10.cc"] <- psi["11.cc"] - psi["10.cc"]
  tau["11.10.ac"] <- psi["11.ac"] - psi["10.ac"]
  tau["11.10.-c"] <- (rho["cc"]/(rho["cc"] + rho["ac"])) * tau["11.10.cc"] +
    (rho["ac"]/(rho["cc"] + rho["ac"])) * tau["11.10.ac"]
  tau["11.00.cc"] <- psi["11.cc"] - psi["00.cc"]
  tau["int.cc"] <- psi["11.cc"] - psi["01.cc"] - psi["10.cc"] + psi["00.cc"]

  Chat <- matrix(0, nrow = length(tau),
                 ncol = length(rho) + length(psi))
  colnames(Chat) <- c(names(rho), names(psi))
  rownames(Chat) <- names(tau)
  Chat["10.00.cc", "10.cc"] <- 1
  Chat["10.00.cc", "00.cc"] <- -1
  Chat["10.00.cn", "10.cn"] <- 1
  Chat["10.00.cn", "00.cn"] <- -1
  Chat["10.00.c-", "10.cc"] <- rho["cc"]/(rho["cc"] + rho["cn"])
  Chat["10.00.c-", "00.cc"] <- -rho["cc"]/(rho["cc"] + rho["cn"])
  Chat["10.00.c-", "10.cn"] <- rho["cc"]/(rho["cc"] + rho["cn"])
  Chat["10.00.c-", "00.cn"] <- -rho["cc"]/(rho["cc"] + rho["cn"])
  Chat["10.00.c-", "cc"] <- tau["10.00.cc"] * (rho["cn"] / (rho["cc"] + rho["cn"]) ^ 2)
  Chat["10.00.c-", "cn"] <- tau["10.00.cn"] * (rho["cc"] / (rho["cc"] + rho["cn"]) ^ 2)
  Chat["01.00.cc", "01.cc"] <- 1
  Chat["01.00.cc", "00.cc"] <- -1
  Chat["01.00.nc", "01.nc"] <- 1
  Chat["01.00.nc", "00.nc"] <- -1
  Chat["01.00.-c", "01.cc"] <- rho["cc"]/(rho["cc"] + rho["nc"])
  Chat["01.00.-c", "00.cc"] <- -rho["cc"]/(rho["cc"] + rho["nc"])
  Chat["01.00.-c", "01.nc"] <- rho["cc"]/(rho["cc"] + rho["nc"])
  Chat["01.00.-c", "00.nc"] <- -rho["cc"]/(rho["cc"] + rho["nc"])
  Chat["01.00.-c", "cc"] <- tau["01.00.cc"] * (rho["nc"] / (rho["cc"] + rho["nc"]) ^ 2)
  Chat["01.00.-c", "nc"] <- tau["01.00.nc"] * (rho["cc"] / (rho["cc"] + rho["nc"]) ^ 2)
  Chat["11.01.cc", "11.cc"] <- 1
  Chat["11.01.cc", "01.cc"] <- -1
  Chat["11.01.ca", "11.ca"] <- 1
  Chat["11.01.ca", "01.ca"] <- -1
  Chat["11.01.c-", "11.cc"] <- rho["cc"]/(rho["cc"] + rho["ca"])
  Chat["11.01.c-", "01.cc"] <- -rho["cc"]/(rho["cc"] + rho["ca"])
  Chat["11.01.c-", "11.ca"] <- rho["cc"]/(rho["cc"] + rho["ca"])
  Chat["11.01.c-", "01.ca"] <- -rho["cc"]/(rho["cc"] + rho["ca"])
  Chat["11.01.c-", "cc"] <- tau["11.01.cc"] * (rho["ca"] / (rho["cc"] + rho["ca"]) ^ 2)
  Chat["11.01.c-", "ca"] <- tau["11.01.ca"] * (rho["cc"] / (rho["cc"] + rho["ca"]) ^ 2)
  Chat["11.10.cc", "11.cc"] <- 1
  Chat["11.10.cc", "10.cc"] <- -1
  Chat["11.10.ac", "11.ac"] <- 1
  Chat["11.10.ac", "10.ac"] <- -1
  Chat["11.10.-c", "11.cc"] <- rho["cc"]/(rho["cc"] + rho["ac"])
  Chat["11.10.-c", "10.cc"] <- -rho["cc"]/(rho["cc"] + rho["ac"])
  Chat["11.10.-c", "11.ac"] <- rho["cc"]/(rho["cc"] + rho["ac"])
  Chat["11.10.-c", "10.ac"] <- -rho["cc"]/(rho["cc"] + rho["ac"])
  Chat["11.10.-c", "cc"] <- tau["11.10.cc"] * (rho["ac"] / (rho["cc"] + rho["ac"]) ^ 2)
  Chat["11.10.-c", "ac"] <- tau["11.10.ac"] * (rho["cc"] / (rho["cc"] + rho["ac"]) ^ 2)
  Chat["11.00.cc", "11.cc"] <- 1
  Chat["11.00.cc", "00.cc"] <- -1
  Chat["int.cc", c("11.cc", "00.cc")] <- 1
  Chat["int.cc", c("10.cc", "01.cc")] <- -1

  return(list(tau = tau, Chat = Chat))
}


iv_int_grad <- function(theta, x) {
  d1 <- x[, 1]
  d2 <- x[, 2]
  z1 <- x[, 3]
  z2 <- x[, 4]
  z11 <- z1 * z2
  z10 <- z1 * (1 - z2)
  z01 <- (1 - z1) * z2
  z00 <- (1 - z1) * (1 - z2)
  d11 <- d1 * d2
  d10 <- d1 * (1 - d2)
  d01 <- (1 - d1) * d2
  d00 <- (1 - d1) * (1 - d2)
  pstr <- expand.grid(c("c", "a", "n"), c("c", "a", "n"))
  pstr_names <- paste0(pstr$Var1, pstr$Var2)
  rho <- c(theta[1:8], 1 - sum(theta[1:8]))
  names(rho) <- pstr_names
  psi <- theta[9:24]
  names(psi) <- c("11.cc", "11.ac", "11.ca", "11.aa",
                  "10.cc", "10.ac", "10.cn", "10.an",
                  "01.cc", "01.nc", "01.na", "01.ca",
                  "00.cc", "00.nc", "00.cn", "00.nn")
  ## "cc" "ac" "nc" "ca" "aa" "na" "cn" "an" "nn"
  rho_mat <- rbind(
    c(1, 1, 0, 1, 1, 0, 0, 0, 0), ## 11.11
    c(0, 0, 1, 0, 0, 1, 0, 0, 0), ## 01.11
    c(0, 0, 0, 0, 0, 0, 1, 1, 0), ## 10.11
    c(0, 1, 0, 0, 1, 0, 0, 0, 0), ## 11.01
    c(1, 0, 1, 1, 0, 1, 0, 0, 0), ## 01.01
    c(0, 0, 0, 0, 0, 0, 0, 1, 0), ## 10.01
    c(0, 0, 0, 1, 1, 0, 0, 0, 0), ## 11.10
    c(0, 0, 0, 0, 0, 1, 0, 0, 0), ## 01.10
    c(1, 1, 0, 0, 0, 0, 1, 1, 0), ## 10.10
    c(0, 0, 0, 0, 1, 0, 0, 0, 0), ## 11.00
    c(0, 0, 0, 1, 0, 1, 0, 0, 0), ## 01.00
    c(0, 1, 0, 0, 0, 0, 0, 1, 0)  ## 10.00
  )
  Zmat <- cbind(z11, z11, z11, z11, z01, z01, z01, z01,
                z10, z10, z10, z10, z00, z00, z00, z00)
  rho_G <- matrix(0, nrow = 28, ncol = 9)
  rho_G[1:12,] <- -colSums(Zmat[, -c(4, 8, 12, 16)]) * rho_mat
  psi_G <- matrix(0, nrow = 28, ncol = 16)
  rownames(psi_G)[13:28] <- c("s11.11", "s01.11", "s10.11", "s00.11",
                              "s11.01", "s01.01", "s10.01", "s00.01",
                              "s11.10", "s01.10", "s10.10", "s00.10",
                              "s11.00", "s01.00", "s10.00", "s00.00")
  rownames(rho_G)[13:28] <- c("s11.11", "s01.11", "s10.11", "s00.11",
                              "s11.01", "s01.01", "s10.01", "s00.01",
                              "s11.10", "s01.10", "s10.10", "s00.10",
                              "s11.00", "s01.00", "s10.00", "s00.00")
  
  colnames(psi_G) <- names(psi)
  colnames(rho_G) <- names(rho)
  psi_mat <- matrix(0, nrow = 16, ncol = 16)

  oneone <- c("11.cc", "11.ca", "11.ac", "11.aa")  
  psi_G["s11.11", oneone] <- -sum(z11) * rho[c("cc", "ca", "ac", "aa")]
  rho_G["s11.11", c("cc", "ca", "ac", "aa")] <- -sum(z11) * psi[oneone]
  psi_G["s01.11", c("01.nc", "01.na")] <- -sum(z11) * rho[c("nc","na")]
  rho_G["s01.11", c("nc","na")] <- -sum(z11) * psi[c("01.nc", "01.na")]
  psi_G["s10.11", c("10.cn", "10.an")] <- -sum(z11) * rho[c("cn","an")]
  rho_G["s10.11", c("cn","an")] <- -sum(z11) * psi[c("10.cn", "10.an")]
  psi_G["s00.11", "00.nn"] <- -sum(z11) * rho["nn"]
  rho_G["s00.11", "nn"] <- -sum(z11) * psi["00.nn"]

  zeroone <- c("01.cc", "01.nc", "01.na", "01.ca")
  psi_G["s11.01", c("11.aa", "11.ac")] <- -sum(z01) * rho[c("aa", "ac")]
  psi_G["s01.01", zeroone] <- -sum(z01) * rho[c("cc", "nc", "na", "ca")]  
  psi_G["s10.01", c("10.an")] <- -sum(z01) * rho["an"]
  psi_G["s00.01", c("00.cn", "00.nn")] <- -sum(z01) * rho[c("cn", "nn")]
  rho_G["s11.01", c("aa", "ac")] <- -sum(z01) * psi[c("11.aa", "11.ac")]
  rho_G["s01.01", c("cc", "nc", "na", "ca")] <- -sum(z01) * psi[zeroone]  
  rho_G["s10.01", "an"] <- -sum(z01) * psi["10.an"]
  rho_G["s00.01", c("cn", "nn")] <- -sum(z01) * psi[c("00.cn", "00.nn")]

  onezero <- c("10.cc", "10.ac", "10.cn", "10.an")
  psi_G["s11.10", c("11.aa", "11.ca")] <- -sum(z10) * rho[c("aa", "ca")]
  psi_G["s01.10", c("01.na")] <- -sum(z10) * rho["na"]
  psi_G["s10.10", onezero] <- -sum(z10) * rho[c("cc", "ac", "cn", "an")]
  psi_G["s00.10", c("00.nn", "00.nc")] <- -sum(z10) * rho[c("nn", "nc")]
  rho_G["s11.10", c("aa", "ca")] <- -sum(z10) * psi[c("11.aa", "11.ca")]
  rho_G["s01.10", c("na")] <- -sum(z10) * psi["01.na"]
  rho_G["s10.10", c("cc", "ac", "cn", "an")] <- -sum(z10) * psi[onezero]
  rho_G["s00.10", c("nn", "nc")] <- -sum(z10) * psi[c("00.nn", "00.nc")]

  zerozero <- c("00.cc", "00.nc", "00.cn", "00.nn")
  psi_G["s11.00", c("11.aa")] <- -sum(z00) * rho["aa"]
  psi_G["s01.00", c("01.na", "01.ca")] <- -sum(z00) * rho[c("na", "ca")]
  psi_G["s10.00", c("10.an", "10.ac")] <- -sum(z00) * rho[c("an", "ac")]
  psi_G["s00.00", zerozero] <- -sum(z00) * rho[c("cc", "nc", "cn", "nn")]
  rho_G["s11.00", c("aa")] <- -sum(z00) * psi["11.aa"]
  rho_G["s01.00", c("na", "ca")] <- -sum(z00) * psi[c("01.na", "01.ca")]
  rho_G["s10.00", c("an", "ac")] <- -sum(z00) * psi[c("10.an", "10.ac")]
  rho_G["s00.00", c("cc", "nc", "cn", "nn")] <- -sum(z00) * psi[zerozero]

  Ghat <- cbind(rho_G[, -9], psi_G)
  return(Ghat)
}

iv_int_Qgrad <- function(theta, x, W) {
  ghat <- colMeans(iv_int_moments(theta = theta, x = x)) 
  Ghat <- iv_int_grad(theta, x)
  Qgrad <- -crossprod(Ghat, W %*% ghat)
  return(Qgrad)
}


iv_init <- function(y, d, z) {
  z11 <- as.logical(z[,1] * z[,2])
  z10 <- as.logical(z[,1] * (1 - z[,2]))
  z01 <- as.logical((1 - z[,1]) * z[,2])
  z00 <- as.logical((1 - z[,1]) * (1 - z[,2]))
  d11 <- d[,1] * d[,2]
  d10 <- d[,1] * (1 - d[,2])
  d01 <- (1 - d[,1]) * d[,2]
  d00 <- (1 - d[,1]) * (1 - d[,2])
  N <- length(y)

  pstr <- expand.grid(c("c", "a", "n"), c("c", "a", "n"))
  pstr_names <- paste0(pstr$Var1, pstr$Var2)
  rho <- rep(NA, times = 9)
  names(rho) <- pstr_names
  rho["cc"] <- mean(d11[z11]) - mean(d11[z10]) - mean(d11[z01]) + mean(d11[z00])
  rho["ca"] <- mean(d11[z10]) - mean(d11[z00])
  rho["cn"] <- mean(d00[z01]) - mean(d00[z11])
  rho["ac"] <- mean(d11[z01]) - mean(d11[z00])
  rho["aa"] <- mean(d11[z00])
  rho["an"] <- mean(d10[z01])
  rho["nc"] <- mean(d00[z10]) - mean(d00[z11])
  rho["na"] <- mean(d01[z10])
  rho["nn"] <- mean(d00[z11])

  psi <- rep(NA, times = 16)
  names(psi) <- c("11.cc", "11.ac", "11.ca", "11.aa",
                  "10.cc", "10.ac", "10.cn", "10.an",
                  "01.cc", "01.nc", "01.na", "01.ca",
                  "00.cc", "00.nc", "00.cn", "00.nn")
  psi["00.nn"] <- mean((y * d00)[z11]) / rho["nn"]
  psi["00.nn"] <- ifelse(rho["nn"] == 0, 0, psi["00.nn"])
  psi["11.aa"] <- mean((y * d11)[z00]) / rho["aa"]
  psi["11.aa"] <- ifelse(rho["aa"] == 0, 0, psi["11.aa"])
  psi["10.an"] <- mean((y * d10)[z01]) / rho["an"]
  psi["10.an"] <- ifelse(rho["an"] == 0, 0, psi["10.an"])
  psi["01.na"] <- mean((y * d01)[z10]) / rho["na"]
  psi["01.na"] <- ifelse(rho["na"] == 0, 0, psi["01.na"])
  psi["01.nc"] <- (mean((y * d01)[z11]) - psi["01.na"] * rho["na"]) / rho["nc"]
  psi["01.nc"] <- ifelse(rho["nc"] == 0, 0, psi["01.nc"])
  psi["10.cn"] <- (mean((y * d10)[z11]) - psi["10.an"] * rho["an"]) / rho["cn"]
  psi["10.cn"] <- ifelse(rho["cn"] == 0, 0, psi["10.cn"])
  psi["11.ac"] <- (mean((y * d11)[z01]) - psi["11.aa"] * rho["aa"]) / rho["ac"]
  psi["11.ac"] <- ifelse(rho["ac"] == 0, 0, psi["11.ac"])
  psi["00.cn"] <- (mean((y * d00)[z01]) - psi["00.nn"] * rho["nn"]) / rho["cn"]
  psi["00.cn"] <- ifelse(rho["cn"] == 0, 0, psi["00.cn"])
  psi["11.ca"] <- (mean((y * d11)[z10]) - psi["11.aa"] * rho["aa"]) / rho["ca"]
  psi["11.ca"] <- ifelse(rho["ca"] == 0, 0, psi["11.ca"])
  psi["00.nc"] <- (mean((y * d00)[z10]) - psi["00.nn"] * rho["nn"]) / rho["nc"]
  psi["00.nc"] <- ifelse(rho["nc"] == 0, 0, psi["00.nc"])
  psi["01.ca"] <- (mean((y * d01)[z00]) - psi["01.na"] * rho["na"]) / rho["ca"]
  psi["01.ca"] <- ifelse(rho["ca"] == 0, 0, psi["01.ca"])
  psi["10.ac"] <- (mean((y * d10)[z00]) - psi["10.an"] * rho["an"]) / rho["ac"]
  psi["10.ac"] <- ifelse(rho["ac"] == 0, 0, psi["10.ac"])
  psi["11.cc"] <- (mean((y * d11)[z11]) - psi["11.ca"] * rho["ca"] -
                     psi["11.ac"] * rho["ac"] -
                     psi["11.aa"] * rho["aa"]) / rho["cc"]
  psi["11.cc"] <- ifelse(rho["cc"] == 0, 0, psi["11.cc"])
  psi["01.cc"] <- (mean((y * d01)[z01]) - psi["01.nc"] * rho["nc"] -
                     psi["01.na"] * rho["na"] -
                     psi["01.ca"] * rho["ca"]) / rho["cc"]
  psi["01.cc"] <- ifelse(rho["cc"] == 0, 0, psi["01.cc"])
  psi["10.cc"] <- (mean((y * d10)[z10]) - psi["10.ac"] * rho["ac"] -
                     psi["10.cn"] * rho["cn"] -
                     psi["10.an"] * rho["an"]) / rho["cc"]
  psi["10.cc"] <- ifelse(rho["cc"] == 0, 0, psi["10.cc"])
  psi["00.cc"] <- (mean((y * d00)[z00]) - psi["00.nc"] * rho["nc"] -
                     psi["00.cn"] * rho["cn"] -
                     psi["00.nn"] * rho["nn"]) / rho["cc"]
  psi["00.cc"] <- ifelse(rho["cc"] == 0, 0, psi["00.cc"])
  
  return(list(rho = rho[1:8], psi = psi))
}
