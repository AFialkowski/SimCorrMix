#' @title Error Loop to Correct Final Correlation of Simulated Variables
#'
#' @description This function attempts to correct the final pairwise correlations of simulated variables to be within \code{epsilon}
#'     of the target correlations.  It updates the intermediate normal correlation iteratively in a loop until either the maximum error
#'     is less than epsilon or the number of iterations exceeds \code{maxit}.  This function would not ordinarily be called directly by
#'     the user.  The function is a modification of  Barbiero & Ferrari's \code{\link[GenOrd]{ordcont}} function in
#'     \code{\link[GenOrd]{GenOrd-package}}.  The \code{\link[GenOrd]{ordcont}} function has been modified in the following ways:
#'
#'     1) It works for continuous, ordinal (r >= 2 categories), and count (regular or zero-inflated, Poisson or Negative Binomial) variables.
#'
#'     2) The initial correlation check has been removed because the intermediate correlation matrix
#'     \code{Sigma} from \code{\link[SimCorrMix]{corrvar}} or \code{\link[SimCorrMix]{corrvar2}} has already been
#'     checked for positive-definiteness and used to generate variables.
#'
#'     3) Eigenvalue decomposition is done on \code{Sigma} to impose the correct intermediate correlations on the normal variables.
#'     If \code{Sigma} is not positive-definite, the negative eigenvalues are replaced with 0.
#'
#'     4) The final positive-definite check has been removed.
#'
#'     5) The intermediate correlation update function was changed to accommodate more situations.
#'
#'     6) Allowing specifications for the sample size and the seed for reproducibility.
#'
#'     The vignette \bold{Variable Types} describes the algorithm used in the error loop.
#'
#' @param n the sample size
#' @param k_cat the number of ordinal (r >= 2 categories) variables
#' @param k_cont the number of continuous variables (these may be regular continuous variables or components of continuous mixture variables)
#' @param k_pois the number of Poisson (regular or zero-inflated) variables
#' @param k_nb the number of Negative Binomial (regular or zero-inflated) variables
#' @param method the method used to generate the continuous variables.  "Fleishman" uses a third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param means a vector of means for the continuous variables
#' @param vars a vector of variances for the continuous variables
#' @param constants a matrix with \code{k_cont} rows, each a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or
#'     c0, c1, c2, c3, c4, c5 (if \code{method} = "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param marginal a list of length equal \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param support a list of length equal \code{k_cat}; the i-th element is a vector of containing the r
#'     ordered support values; if not provided, the default is for the i-th element to be the vector 1, ..., r
#' @param lam a vector of lambda (mean > 0) constants for the Poisson variables (see \code{stats::dpois}); the order should be
#'     1st regular Poisson variables, 2nd zero-inflated Poisson variables
#' @param p_zip a vector of probabilities of structural zeros (not including zeros from the Poisson distribution) for the zero-inflated
#'     Poisson variables (see \code{VGAM::dzipois})
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{stats::dnbinom}); the order should be
#'     1st regular NB variables, 2nd zero-inflated NB variables
#' @param mu a vector of mean parameters for the NB variables; order the same as in \code{size}; for zero-inflated NB this refers to
#'     the mean of the NB distribution (see \code{VGAM::dzinegbin})
#' @param p_zinb a vector of probabilities of structural zeros (not including zeros from the NB distribution) for the zero-inflated NB variables
#'     (see \code{VGAM::dzinegbin})
#' @param seed the seed value for random number generation
#' @param epsilon the maximum acceptable error between the final and target pairwise correlation; smaller epsilons take more time
#' @param maxit the maximum number of iterations to use to find the intermediate correlation; the
#'     correction loop stops when either the iteration number passes \code{maxit} or \code{epsilon} is reached
#' @param rho0 the target correlation matrix
#' @param Sigma the intermediate correlation matrix previously used in \code{\link[SimCorrMix]{corrvar}}
#'     or \code{\link[SimCorrMix]{corrvar2}}
#' @param rho_calc the final correlation matrix calculated in \code{\link[SimCorrMix]{corrvar}}
#'     or \code{\link[SimCorrMix]{corrvar2}} before execution of \code{\link[SimCorrMix]{corr_error}}
#' @importFrom VGAM qzipois qzinegbin
#' @export
#' @keywords error correlation
#' @seealso \code{\link[SimCorrMix]{corrvar}}, \code{\link[SimCorrMix]{corrvar2}}
#'
#' @return A list with the following components:
#' @return \code{Sigma} the intermediate MVN correlation matrix resulting from the error loop
#' @return \code{rho_calc} the calculated final correlation matrix generated from Sigma
#' @return \code{Y_cat} the ordinal variables
#' @return \code{Y} the continuous (mean 0, variance 1) variables
#' @return \code{Y_cont} the continuous variables with desired mean and variance
#' @return \code{Y_pois} the Poisson variables
#' @return \code{Y_nb} the Negative Binomial variables
#' @return \code{niter} a matrix containing the number of iterations required for each variable pair
#' @references Please see references for \code{\link[SimCorrMix]{SimCorrMix}}.
#'
corr_error <- function(n = 10000, k_cat = 0, k_cont = 0, k_pois = 0, k_nb = 0,
                       method = c("Fleishman", "Polynomial"), means = NULL,
                       vars = NULL, constants = NULL, marginal = list(),
                       support = list(), lam = NULL, p_zip = 0,
                       size = NULL, mu = NULL, p_zinb = 0, seed = 1234,
                       epsilon = 0.001, maxit = 1000, rho0 = NULL,
                       Sigma = NULL, rho_calc = NULL) {
  k <- k_cat + k_cont + k_pois + k_nb
  if (k_pois > 0 & length(p_zip) < k_pois)
    p_zip <- c(rep(0, k_pois - length(p_zip)), p_zip)
  if (k_nb > 0 & length(p_zinb) < k_nb)
    p_zinb <- c(rep(0, k_nb - length(p_zinb)), p_zinb)
  niter <- matrix(0, k, k)
  Sigmaold <- Sigma
  Y_cat <- NULL
  Y <- NULL
  Y_cont <- NULL
  Y_pois <- NULL
  Y_nb <- NULL
  for (q in 1:(k - 1)) {
    for (r in (q + 1):k) {
      if (rho0[q, r] == 0) {
        Sigma[q, r] <- 0
      } else {
        it <- 0
        while (abs(rho_calc[q, r] - rho0[q, r]) > epsilon & (it < maxit)) {
          if (rho0[q, r] * (rho0[q, r]/rho_calc[q, r]) <= -1) {
            Sigma[q, r] <- Sigmaold[q, r] *
              (1 + 0.1 * (1 - Sigmaold[q, r]) *
                 -sign(rho0[q, r] - rho_calc[q, r]))
          }
          if (rho0[q, r] * (rho0[q, r]/rho_calc[q, r]) >= 1) {
            Sigma[q, r] <- Sigmaold[q, r] *
              (1 + 0.1 * (1 - Sigmaold[q, r]) *
                 sign(rho0[q, r] - rho_calc[q, r]))
          }
          if ((rho0[q, r] * (rho0[q, r]/rho_calc[q, r]) > -1) &
              (rho0[q, r] * (rho0[q, r]/rho_calc[q, r]) < 1)) {
            Sigma[q, r] <- Sigmaold[q, r] * (rho0[q, r]/rho_calc[q, r])
          }
          if (Sigma[q, r] > 1) Sigma[q, r] <- 1
          if (Sigma[q, r] < -1) Sigma[q, r] <- -1
          Sigma[r, q] <- Sigma[q, r]
          eig <- eigen(Sigma, symmetric = TRUE)
          sqrteigval <- diag(sqrt(pmax(eig$values, 0)))
          eigvec <- eig$vectors
          fry <- eigvec %*% sqrteigval
          set.seed(seed)
          X <- matrix(rnorm(ncol(Sigma) * n), n)
          X <- scale(X, TRUE, FALSE)
          X <- X %*% svd(X, nu = 0)$v
          X <- scale(X, FALSE, TRUE)
          X <- fry %*% t(X)
          X <- t(X)
          if (k_cat > 0) {
            X_cat <- X[, 1:k_cat, drop = FALSE]
            Y_cat <- matrix(1, nrow = n, ncol = ncol(X_cat))
            for (i in 1:ncol(X_cat)) {
              Y_cat[, i] <- as.integer(cut(X_cat[, i],
                                       breaks = c(min(X_cat[, i]) - 1,
                                                  qnorm(marginal[[i]]),
                                                  max(X_cat[, i])  +  1)))
              Y_cat[, i] <- support[[i]][Y_cat[, i]]
            }
          }
          if (k_cont > 0) {
            X_cont <- X[, (k_cat + 1):(k_cat + k_cont), drop = FALSE]
            Y_cont <- matrix(1, nrow = n, ncol = ncol(X_cont))
            Y <- matrix(1, nrow = n, ncol = ncol(X_cont))
            for (i in 1:ncol(X_cont)) {
              if (method == "Fleishman") {
                Y[, i] <- constants[i, 1] +
                  constants[i, 2] * X_cont[, i] +
                  constants[i, 3] * X_cont[, i]^2 +
                  constants[i, 4] * X_cont[, i]^3
              }
              if (method == "Polynomial") {
                Y[, i] <- constants[i, 1] +
                  constants[i, 2] * X_cont[, i] +
                  constants[i, 3] * X_cont[, i]^2 +
                  constants[i, 4] * X_cont[, i]^3 +
                  constants[i, 5] * X_cont[, i]^4 +
                  constants[i, 6] * X_cont[, i]^5
              }
              Y_cont[, i] <- means[i] + sqrt(vars[i]) * Y[, i]
            }
          }
          if (k_pois > 0) {
            X_pois <- X[, (k_cat + k_cont + 1):(k_cat + k_cont + k_pois),
                        drop = FALSE]
            Y_pois <- matrix(1, nrow = n, ncol = ncol(X_pois))
            for (i in 1:ncol(X_pois)) {
              Y_pois[, i] <- qzipois(pnorm(X_pois[, i]), lam[i], p_zip[i])
            }
          }
          if (k_nb > 0) {
            X_nb <- X[, (k_cat + k_cont + k_pois + 1):ncol(X), drop = FALSE]
            Y_nb <- matrix(1, nrow = n, ncol = ncol(X_nb))
            for (i in 1:ncol(X_nb)) {
              Y_nb[, i] <- qzinegbin(pnorm(X_nb[, i]), size = size[i], munb = mu[i],
                                     pstr0 = p_zinb[i])
            }
          }
          rho_calc <- cor(cbind(Y_cat, Y_cont, Y_pois, Y_nb))
          Sigmaold <- Sigma
          it <- it + 1
        }
        niter[q, r] <- it
        niter[r, q] <- it
      }
    }
  }
  rho_calc <- cor(cbind(Y_cat, Y_cont, Y_pois, Y_nb))
  return(list(Sigma = Sigma, rho_calc = rho_calc, Y_cat = Y_cat, Y = Y,
              Y_cont = Y_cont, Y_pois = Y_pois, Y_nb = Y_nb, niter = niter))
}
