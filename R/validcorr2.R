#' @title Determine Correlation Bounds for Ordinal, Continuous, Poisson, and/or Negative Binomial Variables: Correlation Method 2
#'
#' @description This function calculates the lower and upper correlation bounds for the given distributions and
#'     checks if a given target correlation matrix \code{rho} is within the bounds.  It should be used before simulation with
#'     \code{\link[SimCorrMix]{corrvar2}}.  However, even if all pairwise correlations fall within the bounds, it is still possible
#'     that the desired correlation matrix is not feasible.  This is particularly true when ordinal variables (\eqn{r \ge 2} categories) are
#'     generated or negative correlations are desired.  Therefore, this function should be used as a general check to eliminate pairwise correlations that are obviously
#'     not reproducible.  It will help prevent errors when executing the simulation.  The \emph{ordering} of the variables in \code{rho}
#'     must be 1st ordinal, 2nd continuous non-mixture, 3rd components of continuous mixture, 4th regular Poisson, 5th zero-inflated
#'     Poisson, 6th regular NB, and 7th zero-inflated NB.  Note that it is possible for \code{k_cat}, \code{k_cont}, \code{k_mix},
#'     \code{k_pois}, and/or \code{k_nb} to be 0.  The target correlations are specified with respect to the components of the continuous
#'     mixture variables.  There are no parameter input checks in order to decrease simulation time.  All inputs should be checked prior to simulation with
#'     \code{\link[SimCorrMix]{validpar}}.
#'
#'     Please see the \bold{Comparison of Correlation Methods 1 and 2} vignette for the differences between the two correlation methods, and
#'     the \bold{Calculation of Correlation Boundaries} vignette for a detailed explanation of how the correlation boundaries are calculated.
#'
#' @section Reasons for Function Errors:
#'     1) The most likely cause for function errors is that no solutions to \code{\link[SimMultiCorrData]{fleish}} or
#'     \code{\link[SimMultiCorrData]{poly}} converged when using \code{\link[SimMultiCorrData]{find_constants}}.  If this happens,
#'     the function will stop.  It may help to first use \code{\link[SimMultiCorrData]{find_constants}} for each continuous variable to
#'     determine if a sixth cumulant correction value is needed.  If the standardized cumulants are obtained from \code{calc_theory},
#'     the user may need to use rounded values as inputs (i.e. \code{skews = round(skews, 8)}).  For example, in order to ensure that skew
#'     is exactly 0 for symmetric distributions.
#'
#'     2) The kurtosis may be outside the region of possible values.  There is an associated lower boundary for kurtosis associated
#'     with a given skew (for Fleishman's method) or skew and fifth and sixth cumulants (for Headrick's method).  Use
#'     \code{\link[SimMultiCorrData]{calc_lower_skurt}} to determine the boundary for a given set of cumulants.
#'
#' @param n the sample size (i.e. the length of each simulated variable; default = 10000)
#' @param k_cat the number of ordinal (r >= 2 categories) variables (default = 0)
#' @param k_cont the number of continuous non-mixture variables (default = 0)
#' @param k_mix the number of continuous mixture variables (default = 0)
#' @param k_pois the number of regular Poisson and zero-inflated Poisson variables (default = 0)
#' @param k_nb the number of regular Negative Binomial and zero-inflated Negative Binomial variables (default = 0)
#' @param method the method used to generate the k_cont non-mixture and k_mix mixture continuous variables.  "Fleishman" uses
#'     Fleishman's third-order polynomial transformation and "Polynomial" uses Headrick's fifth-order transformation.
#' @param means a vector of means for the k_cont non-mixture and k_mix mixture continuous variables
#'     (i.e. \code{rep(0, (k_cont + k_mix))})
#' @param vars a vector of variances for the k_cont non-mixture and k_mix mixture continuous variables
#'     (i.e. \code{rep(1, (k_cont + k_mix))})
#' @param skews a vector of skewness values for the \code{k_cont} non-mixture continuous variables
#' @param skurts a vector of standardized kurtoses (kurtosis - 3, so that normal variables have a value of 0)
#'     for the \code{k_cont} non-mixture continuous variables
#' @param fifths a vector of standardized fifth cumulants for the \code{k_cont} non-mixture continuous variables
#'     (not necessary for \code{method} = "Fleishman")
#' @param sixths a vector of standardized sixth cumulants for the \code{k_cont} non-mixture continuous variables
#'     (not necessary for \code{method} = "Fleishman")
#' @param Six a list of vectors of sixth cumulant correction values for the \code{k_cont} non-mixture continuous variables
#'     if no valid PDF constants are found, \cr ex: \code{Six = list(seq(0.01, 2, 0.01), seq(1, 10, 0.5))};
#'     if no correction is desired for variable \eqn{Y_{cont_i}}, set set the i-th list component equal to \code{NULL};
#'     if no correction is desired for any of the \eqn{Y_{cont}} keep as \code{Six = list()}
#'     (not necessary for \code{method} = "Fleishman")
#' @param mix_pis a list of length \code{k_mix} with i-th component a vector of mixing probabilities that sum to 1 for component distributions of \eqn{Y_{mix_i}}
#' @param mix_mus a list of length \code{k_mix} with i-th component a vector of means for component distributions of \eqn{Y_{mix_i}}
#' @param mix_sigmas a list of length \code{k_mix} with i-th component a vector of standard deviations for component distributions of \eqn{Y_{mix_i}}
#' @param mix_skews a list of length \code{k_mix} with i-th component a vector of skew values for component distributions of \eqn{Y_{mix_i}}
#' @param mix_skurts a list of length \code{k_mix} with i-th component a vector of standardized kurtoses for component distributions of \eqn{Y_{mix_i}}
#' @param mix_fifths a list of length \code{k_mix} with i-th component a vector of standardized fifth cumulants for component distributions of \eqn{Y_{mix_i}}
#'     (not necessary for \code{method} = "Fleishman")
#' @param mix_sixths a list of length \code{k_mix} with i-th component a vector of standardized sixth cumulants for component distributions of \eqn{Y_{mix_i}}
#'     (not necessary for \code{method} = "Fleishman")
#' @param mix_Six a list of length \code{k_mix} with i-th component a list of vectors of sixth cumulant correction values
#'     for component distributions of \eqn{Y_{mix_i}}; use \code{NULL} if no correction is desired for a given component or
#'     mixture variable; if no correction is desired for any of the \eqn{Y_{mix}} keep as \code{mix_Six = list()}
#'     (not necessary for \code{method} = "Fleishman")
#' @param marginal a list of length equal to \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1);
#'     for binary variables, these should be input the same as for ordinal variables with more than 2 categories (i.e. the user-specified
#'     probability is the probability of the 1st category, which has the smaller support value)
#' @param lam a vector of lambda (> 0) constants for the Poisson variables (see \code{\link[stats]{dpois}}); the order should be
#'     1st regular Poisson variables, 2nd zero-inflated Poisson variables
#' @param p_zip a vector of probabilities of structural zeros (not including zeros from the Poisson distribution) for the
#'     zero-inflated Poisson variables (see \code{\link[VGAM]{dzipois}}); if \code{p_zip} = 0, \eqn{Y_{pois}} has a regular Poisson
#'     distribution; if \code{p_zip} is in (0, 1), \eqn{Y_{pois}} has a zero-inflated Poisson distribution;
#'     if \code{p_zip} is in \code{(-(exp(lam) - 1)^(-1), 0)}, \eqn{Y_{pois}} has a zero-deflated Poisson distribution and \code{p_zip}
#'     is not a probability; if \code{p_zip = -(exp(lam) - 1)^(-1)}, \eqn{Y_{pois}} has a positive-Poisson distribution
#'     (see \code{\link[VGAM]{dpospois}}); if \code{length(p_zip) < length(lam)}, the missing values are set to 0 (and ordered 1st)
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{dnbinom}}); the order should be
#'     1st regular NB variables, 2nd zero-inflated NB variables
#' @param prob a vector of success probability parameters for the NB variables; order the same as in \code{size}
#' @param mu a vector of mean parameters for the NB variables (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL); order the same as in \code{size}
#' @param p_zinb a vector of probabilities of structural zeros (not including zeros from the NB distribution) for the zero-inflated NB variables
#'     (see \code{\link[VGAM]{dzinegbin}}); if \code{p_zinb} = 0, \eqn{Y_{nb}} has a regular NB distribution;
#'     if \code{p_zinb} is in \code{(-prob^size/(1 - prob^size),} \code{0)}, \eqn{Y_{nb}} has a zero-deflated NB distribution and \code{p_zinb}
#'     is not a probability; if \code{p_zinb = -prob^size/(1 - prob^size)}, \eqn{Y_{nb}} has a positive-NB distribution (see
#'     \code{\link[VGAM]{dposnegbin}}); if \code{length(p_zinb) < length(size)}, the missing values are set to 0 (and ordered 1st)
#' @param pois_eps a vector of length \code{k_pois} containing total cumulative probability truncation values; if none are provided,
#'     the default is 0.0001 for each variable
#' @param nb_eps a vector of length \code{k_nb} containing total cumulative probability truncation values; if none are provided,
#'     the default is 0.0001 for each variable
#' @param rho the target correlation matrix which must be ordered
#'     \emph{1st ordinal, 2nd continuous non-mixture, 3rd components of continuous mixtures, 4th regular Poisson, 5th zero-inflated Poisson,
#'     6th regular NB, 7th zero-inflated NB}; note that \code{rho} is specified in terms of the components of \eqn{Y_{mix}}
#' @param seed the seed value for random number generation (default = 1234)
#' @param use.nearPD TRUE to convert \code{rho} to the nearest positive definite matrix with \code{Matrix::nearPD} if necessary
#' @import SimMultiCorrData
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @import nleqslv
#' @import BB
#' @importFrom Matrix nearPD
#' @importFrom VGAM dzipois dzinegbin
#' @export
#' @keywords correlation bounds method2
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimCorrMix]{corrvar2}}, \code{\link[SimCorrMix]{validpar}}
#' @return A list with components:
#' @return \code{rho} the target correlation matrix, which will differ from the supplied matrix (if provided) if it was converted to
#'     the nearest positive-definite matrix
#' @return \code{L_rho} the lower correlation bound
#' @return \code{U_rho} the upper correlation bound
#' @return If continuous variables are desired, additional components are:
#' @return \code{constants} the calculated constants
#' @return \code{sixth_correction} a vector of the sixth cumulant correction values
#' @return \code{valid.pdf} a vector with i-th component equal to "TRUE" if variable Y_i has a valid power method PDF, else "FALSE"
#' @return If a target correlation matrix \code{rho} is provided, each pairwise correlation is checked to see if it is within the lower and upper
#' bounds.  If the correlation is outside the bounds, the indices of the variable pair are given.
#' @return \code{valid.rho} TRUE if all entries of \code{rho} are within the bounds, else FALSE
#' @references Please see \code{\link[SimCorrMix]{corrvar2}} and \code{\link[SimCorrMix]{validcorr}} for references.
#'
#' @examples \dontrun{
#'
#' # 2 continuous mixture, 1 binary, 1 zero-inflated Poisson, and
#' # 1 zero-inflated NB variable
#' n <- 10000
#' seed <- 1234
#'
#' # Mixture variables: Normal mixture with 2 components;
#' # mixture of Logistic(0, 1), Chisq(4), Beta(4, 1.5)
#' # Find cumulants of components of 2nd mixture variable
#' L <- calc_theory("Logistic", c(0, 1))
#' C <- calc_theory("Chisq", 4)
#' B <- calc_theory("Beta", c(4, 1.5))
#'
#' skews <- skurts <- fifths <- sixths <- NULL
#' Six <- list()
#' mix_pis <- list(c(0.4, 0.6), c(0.3, 0.2, 0.5))
#' mix_mus <- list(c(-2, 2), c(L[1], C[1], B[1]))
#' mix_sigmas <- list(c(1, 1), c(L[2], C[2], B[2]))
#' mix_skews <- list(rep(0, 2), c(L[3], C[3], B[3]))
#' mix_skurts <- list(rep(0, 2), c(L[4], C[4], B[4]))
#' mix_fifths <- list(rep(0, 2), c(L[5], C[5], B[5]))
#' mix_sixths <- list(rep(0, 2), c(L[6], C[6], B[6]))
#' mix_Six <- list(list(NULL, NULL), list(1.75, NULL, 0.03))
#' Nstcum <- calc_mixmoments(mix_pis[[1]], mix_mus[[1]], mix_sigmas[[1]],
#'   mix_skews[[1]], mix_skurts[[1]], mix_fifths[[1]], mix_sixths[[1]])
#' Mstcum <- calc_mixmoments(mix_pis[[2]], mix_mus[[2]], mix_sigmas[[2]],
#'   mix_skews[[2]], mix_skurts[[2]], mix_fifths[[2]], mix_sixths[[2]])
#' means <- c(Nstcum[1], Mstcum[1])
#' vars <- c(Nstcum[2]^2, Mstcum[2]^2)
#'
#' marginal <- list(0.3)
#' support <- list(c(0, 1))
#' lam <- 0.5
#' p_zip <- 0.1
#' pois_eps <- 0.0001
#' size <- 2
#' prob <- 0.75
#' p_zinb <- 0.2
#' nb_eps <- 0.0001
#'
#' k_cat <- k_pois <- k_nb <- 1
#' k_cont <- 0
#' k_mix <- 2
#' Rey <- matrix(0.39, 8, 8)
#' diag(Rey) <- 1
#' rownames(Rey) <- colnames(Rey) <- c("O1", "M1_1", "M1_2", "M2_1", "M2_2",
#'   "M2_3", "P1", "NB1")
#'
#' # set correlation between components of the same mixture variable to 0
#' Rey["M1_1", "M1_2"] <- Rey["M1_2", "M1_1"] <- 0
#' Rey["M2_1", "M2_2"] <- Rey["M2_2", "M2_1"] <- Rey["M2_1", "M2_3"] <- 0
#' Rey["M2_3", "M2_1"] <- Rey["M2_2", "M2_3"] <- Rey["M2_3", "M2_2"] <- 0
#'
#' # check parameter inputs
#' validpar(k_cat, k_cont, k_mix, k_pois, k_nb, "Polynomial", means,
#'   vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
#'   mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support,
#'   lam, p_zip, size, prob, mu = NULL, p_zinb, pois_eps, nb_eps, Rey)
#'
#' # check to make sure Rey is within the feasible correlation boundaries
#' validcorr2(n, k_cat, k_cont, k_mix, k_pois, k_nb, "Polynomial", means,
#'   vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
#'   mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
#'   lam, p_zip, size, prob, mu = NULL, p_zinb, pois_eps, nb_eps, Rey, seed)
#' }
validcorr2 <- function(n = 10000, k_cat = 0, k_cont = 0, k_mix = 0, k_pois = 0,
                       k_nb = 0, method = c("Fleishman", "Polynomial"),
                       means =  NULL, vars =  NULL, skews =  NULL,
                       skurts =  NULL, fifths =  NULL, sixths =  NULL,
                       Six = list(), mix_pis = list(), mix_mus = list(),
                       mix_sigmas = list(), mix_skews =  list(),
                       mix_skurts =  list(), mix_fifths =  list(),
                       mix_sixths =  list(), mix_Six = list(),
                       marginal = list(), lam  =  NULL, p_zip = 0, size = NULL,
                       prob = NULL, mu = NULL, p_zinb = 0, pois_eps = 0.0001,
                       nb_eps = 0.0001, rho = NULL, seed = 1234,
                       use.nearPD = TRUE) {
  if (k_pois > 0) {
    if (length(p_zip) < k_pois)
      p_zip <- c(rep(0, k_pois - length(p_zip)), p_zip)
    if (length(pois_eps) < k_pois)
      pois_eps <- rep(0.0001, k_pois)
  }
  if (k_nb > 0) {
    if (length(prob) > 0)
      mu <- size * (1 - prob)/prob
    if (length(p_zinb) < k_nb)
      p_zinb <- c(rep(0, k_nb - length(p_zinb)), p_zinb)
    if (length(nb_eps) < k_nb)
      nb_eps <- rep(0.0001, k_nb)
  }
  if (!is.null(rho)) {
    if (!isSymmetric(rho) | !all(diag(rho) == 1))
      stop("Correlation matrix not valid! Check symmetry and diagonal values.")
    if (min(eigen(rho, symmetric = TRUE)$values) < 0) {
      if (use.nearPD == TRUE) {
        message("Target correlation matrix is not positive definite.
                Nearest positive definite matrix is used!")
        rho <- as.matrix(nearPD(rho, corr = T, keepDiag = T)$mat)
      } else {
        stop("Target correlation matrix is not positive definite.
             Set use.nearPD = TRUE to use nearest positive definite matrix.")
      }
    }
  }
  csame.dist <- NULL
  msame.dist <- NULL
  if (length(skews) >= 2) {
    for (i in 2:length(skews)) {
      if (skews[i] %in% skews[1:(i - 1)]) {
        csame <- which(skews[1:(i - 1)] == skews[i])
        for (j in 1:length(csame)) {
          if (method == "Polynomial") {
            if ((skurts[i] == skurts[csame[j]]) &
                (fifths[i] == fifths[csame[j]]) &
                (sixths[i] == sixths[csame[j]])) {
              csame.dist <- rbind(csame.dist, c(csame[j], i))
              break
            }
          }
          if (method == "Fleishman") {
            if (skurts[i] == skurts[csame[j]]) {
              csame.dist <- rbind(csame.dist, c(csame[j], i))
              break
            }
          }
        }
      }
    }
  }
  mix_skews2 <- NULL
  mix_skurts2 <- NULL
  mix_fifths2 <- NULL
  mix_sixths2 <- NULL
  mix_Six2 <- list()
  if (length(mix_pis) >= 1) {
    k.comp <- c(0, cumsum(unlist(lapply(mix_pis, length))))
    mix_pis2 <- unlist(mix_pis)
    mix_mus2 <- unlist(mix_mus)
    mix_sigmas2 <- unlist(mix_sigmas)
    mix_skews2 <- unlist(mix_skews)
    mix_skurts2 <- unlist(mix_skurts)
    if (method == "Polynomial") {
      mix_fifths2 <- unlist(mix_fifths)
      mix_sixths2 <- unlist(mix_sixths)
      if (length(mix_Six) > 0) {
        if (class(mix_Six[[1]]) == "numeric") mix_Six2 <- mix_Six
        if (class(mix_Six[[1]]) == "list") mix_Six2 <- do.call(append, mix_Six)
      }
    }
    for (i in 1:length(mix_skews2)) {
      msame.dist2 <- NULL
      if (length(skews) >= 1) {
        if (mix_skews2[i] %in% skews) {
          msame <- which(skews == mix_skews2[i])
          for (j in 1:length(msame)) {
            if (method == "Polynomial") {
              if ((mix_skurts2[i] == skurts[msame[j]]) &
                  (mix_fifths2[i] == fifths[msame[j]]) &
                  (mix_sixths2[i] == sixths[msame[j]])) {
                msame.dist2 <- c(msame[j], i)
                break
              }
            }
            if (method == "Fleishman") {
              if (mix_skurts2[i] == skurts[msame[j]]) {
                msame.dist2 <- c(msame[j], i)
                break
              }
            }
          }
        }
      }
      if (is.null(msame.dist2) & i >= 2) {
        if (mix_skews2[i] %in% mix_skews2[1:(i-1)]) {
          msame <- which(mix_skews2[1:(i - 1)] == mix_skews2[i])
          for (j in 1:length(msame)) {
            if (method == "Polynomial") {
              if ((mix_skurts2[i] == mix_skurts2[msame[j]]) &
                  (mix_fifths2[i] == mix_fifths2[msame[j]]) &
                  (mix_sixths2[i] == mix_sixths2[msame[j]])) {
                msame.dist2 <- c(k_cont + msame[j], i)
                break
              }
            }
            if (method == "Fleishman") {
              if (mix_skurts2[i] == mix_skurts2[msame[j]]) {
                msame.dist2 <- c(k_cont + msame[j], i)
                break
              }
            }
          }
        }
      }
      msame.dist <- rbind(msame.dist, msame.dist2)
    }
  }
  if ((k_cont + k_mix) >= 1) {
    SixCorr <- numeric(k_cont + length(mix_skews2))
    Valid.PDF <- numeric(k_cont + length(mix_skews2))
    if (method == "Fleishman") {
      constants <- matrix(NA, nrow = (k_cont + length(mix_skews2)), ncol = 4)
      colnames(constants) <- c("c0", "c1", "c2", "c3")
    }
    if (method == "Polynomial") {
      constants <- matrix(NA, nrow = (k_cont + length(mix_skews2)), ncol = 6)
      colnames(constants) <- c("c0", "c1", "c2", "c3", "c4", "c5")
    }
    set.seed(seed)
    X_cont <- matrix(rnorm((k_cont + length(mix_skews2)) * n), n,
                     (k_cont + length(mix_skews2)))
    X_cont <- scale(X_cont, TRUE, FALSE)
    X_cont <- X_cont %*% svd(X_cont, nu = 0)$v
    X_cont <- scale(X_cont, FALSE, TRUE)
  }
  if (k_cont >= 1) {
    for (i in 1:k_cont) {
      if (!is.null(csame.dist)) {
        rind <- which(csame.dist[, 2] == i)
        if (length(rind) > 0) {
          constants[i, ] <- constants[csame.dist[rind, 1], ]
          SixCorr[i] <- SixCorr[csame.dist[rind, 1]]
          Valid.PDF[i] <- Valid.PDF[csame.dist[rind, 1]]
        }
      }
      if (sum(is.na(constants[i, ])) > 0) {
        if (length(Six) == 0) Six2 <- NULL else
          Six2 <- Six[[i]]
        cons <-
          suppressWarnings(find_constants(method = method, skews = skews[i],
            skurts = skurts[i], fifths = fifths[i], sixths = sixths[i],
            Six = Six2, n = 25, seed = seed))
        if (length(cons) == 1 | is.null(cons)) {
          stop(paste("Constants can not be found for continuous variable ", i,
                     ".", sep = ""))
        }
        con_solution <- cons$constants
        SixCorr[i] <- ifelse(is.null(cons$SixCorr1), NA, cons$SixCorr1)
        Valid.PDF[i] <- cons$valid
        constants[i, ] <- con_solution
      }
    }
  }
  if (k_mix >= 1) {
    for (i in 1:length(mix_skews2)) {
      if (!is.null(msame.dist)) {
        rind <- which(msame.dist[, 2] == i)
        if (length(rind) > 0) {
          constants[(k_cont + i), ] <- constants[msame.dist[rind, 1], ]
          SixCorr[k_cont + i] <- SixCorr[msame.dist[rind, 1]]
          Valid.PDF[k_cont + i] <- Valid.PDF[msame.dist[rind, 1]]
        }
      }
      if (sum(is.na(constants[(k_cont + i), ])) > 0) {
        if (length(mix_Six2) == 0) Six2 <- NULL else
          Six2 <- mix_Six2[[i]]
        cons <-
          suppressWarnings(find_constants(method = method,
            skews = mix_skews2[i], skurts = mix_skurts2[i],
            fifths = mix_fifths2[i], sixths = mix_sixths2[i], Six = Six2,
            n = 25, seed = seed))
        if (length(cons) == 1 | is.null(cons)) {
          stop(paste("Constants can not be found for component variable ", i,
                     ".", sep = ""))
        }
        con_solution <- cons$constants
        SixCorr[k_cont + i] <- ifelse(is.null(cons$SixCorr1), NA,
                                      cons$SixCorr1)
        Valid.PDF[k_cont + i] <- cons$valid
        constants[(k_cont + i), ] <- con_solution
      }
    }
  }
  if ((k_cont + length(mix_skews2)) > 0) {
    Y <- matrix(1, nrow = n, ncol = k_cont + length(mix_skews2))
    for (i in 1:(k_cont + length(mix_skews2))) {
      if (method == "Fleishman") {
        Y[, i] <- constants[i, 1] + constants[i, 2] * X_cont[, i] +
          constants[i, 3] * X_cont[, i]^2 + constants[i, 4] * X_cont[, i]^3
      }
      if (method == "Polynomial") {
        Y[, i] <- constants[i, 1] + constants[i, 2] * X_cont[, i] +
          constants[i, 3] * X_cont[, i]^2 + constants[i, 4] * X_cont[, i]^3 +
          constants[i, 5] * X_cont[, i]^4 + constants[i, 6] * X_cont[, i]^5
      }
    }
    Y_cont <- NULL
    Y_mix <- NULL
    if (k_cont > 0) {
      Y_cont <- matrix(1, n, k_cont)
      means2 <- means[1:k_cont]
      vars2 <- vars[1:k_cont]
      for (i in 1:k_cont) {
        Y_cont[, i] <- means2[i] + sqrt(vars2[i]) * Y[, i]
      }
    }
    if (length(mix_skews2) > 0) {
      Y_mix <- matrix(1, n, length(mix_skews2))
      for (i in 1:length(mix_skews2)) {
        Y_mix[, i] <- mix_mus2[i] + mix_sigmas2[i] * Y[, (k_cont + i)]
      }
    }
    Y <- cbind(Y_cont, Y_mix)
  }
  k_cont <- k_cont + length(mix_skews2)
  k <- k_cat + k_cont + k_pois + k_nb
  L_sigma <- diag(k)
  U_sigma <- diag(k)
  set.seed(seed)
  u <- runif(n, 0, 1)
  set.seed(seed + 1)
  rnorms <- matrix(rnorm(2 * n, 0, 1), ncol = 2)
  rnorms <- scale(rnorms, TRUE, FALSE)
  rnorms <- rnorms %*% svd(rnorms, nu = 0)$v
  rnorms <- scale(rnorms, FALSE, TRUE)
  if (k_pois > 0 | k_nb > 0) {
    max_support <- maxcount_support(k_pois = k_pois, k_nb = k_nb, lam = lam,
      pois_eps = pois_eps, p_zip = p_zip, size = size, mu = mu,
      nb_eps = nb_eps, p_zinb = p_zinb)
  }
  if (k_pois > 0) {
    pois_max <- max_support[max_support$Distribution == "Poisson", ]
    pois_support <- list()
    pois_prob <- list()
    pois_marg <- list()
    for (i in 1:k_pois) {
      pois_support[[i]] <- 0:pois_max[i, 3]
      pois_prob[[i]] <- dzipois(0:pois_max[i, 3], lam[i], p_zip[i])
      pois_prob[[i]][pois_max[i, 3] + 1] <-
        1 - sum(pois_prob[[i]][1:pois_max[i, 3]])
      pois_marg[[i]] <- cumsum(pois_prob[[i]])
      pois_marg[[i]] <- pois_marg[[i]][-(pois_max[i, 3] + 1)]
    }
  }
  if (k_nb > 0) {
    nb_max <- max_support[max_support$Distribution == "Neg_Bin", ]
    nb_support <- list()
    nb_prob <- list()
    nb_marg <- list()
    for (i in 1:k_nb) {
      nb_support[[i]] <- 0:nb_max[i, 3]
      nb_prob[[i]] <- dzinegbin(0:nb_max[i, 3], size = size[i], munb = mu[i],
                                pstr0 = p_zinb[i])
      nb_prob[[i]][nb_max[i, 3] + 1] <- 1 - sum(nb_prob[[i]][1:nb_max[i, 3]])
      nb_marg[[i]] <- cumsum(nb_prob[[i]])
      nb_marg[[i]] <- nb_marg[[i]][-(nb_max[i, 3] + 1)]
    }
  }
  for (q in 1:(k - 1)) {
    for (r in (q + 1):k) {
      if (q >= 1 & q <= k_cat & r >= 1 & r <= k_cat) {
        marg1 <- marginal[[q]]
        marg2 <- marginal[[r]]
        if (length(marg1) == 1 & length(marg2) == 1) {
          p1 <- marg1
          q1 <- 1 - marg1
          p2 <- marg2
          q2 <- 1 - marg2
          L_sigma[q, r] <- L_sigma[r, q] <- max(-sqrt((p1 * p2)/(q1 * q2)),
                                                -sqrt((q1 * q2)/(p1 * p2)))
          U_sigma[q, r] <- U_sigma[r, q] <- min(sqrt((p1 * q2)/(q1 * p2)),
                                                sqrt((q1 * p2)/(p1 * q2)))
        } else {
          n1 <- rnorms[, 1]
          n2 <- rnorms[, 2]
          nord1 <- numeric(length(n1))
          nord2 <- numeric(length(n2))
          for (i in 1:length(marg1)) {
            if (i != length(marg1)) {
              q1 <- qnorm(marg1[i])
              q2 <- qnorm(marg1[i + 1])
              nord1[(q1 < n1) & (n1 <= q2)] <- i
            } else {
              nord1[n1 > qnorm(marg1[i])] <- i
            }
          }
          for (i in 1:length(marg2)) {
            if (i != length(marg2)) {
              q1 <- qnorm(marg2[i])
              q2 <- qnorm(marg2[i + 1])
              nord2[(q1 < n2) & (n2 <= q2)] <- i
            } else {
              nord2[n2 > qnorm(marg2[i])] <- i
            }
          }
          nord1 <- nord1 + 1
          nord2 <- nord2 + 1
          L_sigma[q, r] <- L_sigma[r, q] <-
            cor(nord1[order(nord1, decreasing = TRUE)], nord2[order(nord2)])
          U_sigma[q, r] <- U_sigma[r, q] <- cor(nord1[order(nord1)],
                                                nord2[order(nord2)])
        }
      }
      if (q >= 1 & q <= k_cat & r >= (k_cat + 1) & r <= (k_cat + k_cont)) {
        n1 <- rnorms[, 1]
        n2 <- Y[, r - k_cat]
        marg1 <- marginal[[q]]
        nord1 <- numeric(length(n1))
        for (i in 1:length(marg1)) {
          if (i != length(marg1)) {
            q1 <- qnorm(marg1[i])
            q2 <- qnorm(marg1[i + 1])
            nord1[(q1 < n1) & (n1 <= q2)] <- i
          } else {
            nord1[n1 > qnorm(marg1[i])] <- i
          }
        }
        L_sigma[q, r] <- L_sigma[r, q] <-
          cor(nord1[order(nord1, decreasing = TRUE)], n2[order(n2)])
        U_sigma[q, r] <- U_sigma[r, q] <- cor(nord1[order(nord1)],
                                              n2[order(n2)])
      }
      if (q >= 1 & q <= k_cat & r >= (k_cat + k_cont + 1) &
          r <= (k_cat + k_cont + k_pois)) {
        marg1 <- marginal[[q]]
        marg2 <- pois_marg[[r - (k_cat + k_cont)]]
        n1 <- rnorms[, 1]
        n2 <- rnorms[, 2]
        nord1 <- numeric(length(n1))
        nord2 <- numeric(length(n2))
        for (i in 1:length(marg1)) {
          if (i != length(marg1)) {
            q1 <- qnorm(marg1[i])
            q2 <- qnorm(marg1[i + 1])
            nord1[(q1 < n1) & (n1 <= q2)] <- i
          } else {
            nord1[n1 > qnorm(marg1[i])] <- i
          }
        }
        for (i in 1:length(marg2)) {
          if (i != length(marg2)) {
            q1 <- qnorm(marg2[i])
            q2 <- qnorm(marg2[i + 1])
            nord2[(q1 < n2) & (n2 <= q2)] <- i
          } else {
            nord2[n2 > qnorm(marg2[i])] <- i
          }
        }
        nord1 <- nord1 + 1
        nord2 <- nord2 + 1
        L_sigma[q, r] <- L_sigma[r, q] <-
          cor(nord1[order(nord1, decreasing = TRUE)], nord2[order(nord2)])
        U_sigma[q, r] <- U_sigma[r, q] <- cor(nord1[order(nord1)],
                                              nord2[order(nord2)])
      }
      if (q >= 1 & q <= k_cat & r >= (k_cat + k_cont + k_pois + 1) &
          r <= (k_cat + k_cont + k_pois + k_nb)) {
        marg1 <- marginal[[q]]
        marg2 <- nb_marg[[r - (k_cat + k_cont + k_pois)]]
        n1 <- rnorms[, 1]
        n2 <- rnorms[, 2]
        nord1 <- numeric(length(n1))
        nord2 <- numeric(length(n2))
        for (i in 1:length(marg1)) {
          if (i != length(marg1)) {
            q1 <- qnorm(marg1[i])
            q2 <- qnorm(marg1[i + 1])
            nord1[(q1 < n1) & (n1 <= q2)] <- i
          } else {
            nord1[n1 > qnorm(marg1[i])] <- i
          }
        }
        for (i in 1:length(marg2)) {
          if (i != length(marg2)) {
            q1 <- qnorm(marg2[i])
            q2 <- qnorm(marg2[i + 1])
            nord2[(q1 < n2) & (n2 <= q2)] <- i
          } else {
            nord2[n2 > qnorm(marg2[i])] <- i
          }
        }
        nord1 <- nord1 + 1
        nord2 <- nord2 + 1
        L_sigma[q, r] <- L_sigma[r, q] <-
          cor(nord1[order(nord1, decreasing = TRUE)], nord2[order(nord2)])
        U_sigma[q, r] <- U_sigma[r, q] <- cor(nord1[order(nord1)],
                                              nord2[order(nord2)])
      }
      if (q >= (k_cat + 1) & q <= (k_cat + k_cont) & r >= (k_cat + 1) &
          r <= (k_cat + k_cont)) {
        n1 <- Y[, q - k_cat]
        n2 <- Y[, r - k_cat]
        L_sigma[q, r] <- L_sigma[r, q] <- cor(n1[order(n1, decreasing = TRUE)],
                                              n2[order(n2)])
        U_sigma[q, r] <- U_sigma[r, q] <- cor(n1[order(n1)], n2[order(n2)])
      }
      if (q >= (k_cat + 1) & q <= (k_cat + k_cont) &
          r >= (k_cat + k_cont + 1) & r <= (k_cat + k_cont + k_pois)) {
        n1 <- Y[, q - k_cat]
        n2 <- rnorms[, 1]
        marg1 <- pois_marg[[r - (k_cat + k_cont)]]
        nord1 <- numeric(length(n2))
        for (i in 1:length(marg1)) {
          if (i != length(marg1)) {
            q1 <- qnorm(marg1[i])
            q2 <- qnorm(marg1[i + 1])
            nord1[(q1 < n2) & (n2 <= q2)] <- i
          } else {
            nord1[n2 > qnorm(marg1[i])] <- i
          }
        }
        L_sigma[q, r] <- L_sigma[r, q] <-
          cor(n1[order(n1, decreasing = TRUE)], nord1[order(nord1)])
        U_sigma[q, r] <- U_sigma[r, q] <- cor(n1[order(n1)],
                                              nord1[order(nord1)])
      }
      if (q >= (k_cat + 1) & q <= (k_cat + k_cont) &
         r >= (k_cat + k_cont + k_pois + 1) &
         r <= (k_cat + k_cont + k_pois + k_nb)) {
        n1 <- Y[, q - k_cat]
        n2 <- rnorms[, 1]
        marg1 <- nb_marg[[r - (k_cat + k_cont + k_pois)]]
        nord1 <- numeric(length(n2))
        for (i in 1:length(marg1)) {
          if (i != length(marg1)) {
            q1 <- qnorm(marg1[i])
            q2 <- qnorm(marg1[i + 1])
            nord1[(q1 < n2) & (n2 <= q2)] <- i
          } else {
            nord1[n2 > qnorm(marg1[i])] <- i
          }
        }
        L_sigma[q, r] <- L_sigma[r, q] <- cor(n1[order(n1, decreasing = TRUE)],
                                              nord1[order(nord1)])
        U_sigma[q, r] <- U_sigma[r, q] <- cor(n1[order(n1)],
                                              nord1[order(nord1)])
      }
      if (q >= (k_cat + k_cont + 1) & q <= (k_cat + k_cont + k_pois) &
          r >= (k_cat + k_cont + 1) & r <= (k_cat + k_cont + k_pois)) {
        marg1 <- pois_marg[[q - (k_cat + k_cont)]]
        marg2 <- pois_marg[[r - (k_cat + k_cont)]]
        n1 <- rnorms[, 1]
        n2 <- rnorms[, 2]
        nord1 <- numeric(length(n1))
        nord2 <- numeric(length(n2))
        for (i in 1:length(marg1)) {
          if (i != length(marg1)) {
            q1 <- qnorm(marg1[i])
            q2 <- qnorm(marg1[i + 1])
            nord1[(q1 < n1) & (n1 <= q2)] <- i
          } else {
            nord1[n1 > qnorm(marg1[i])] <- i
          }
        }
        for (i in 1:length(marg2)) {
          if (i != length(marg2)) {
            q1 <- qnorm(marg2[i])
            q2 <- qnorm(marg2[i + 1])
            nord2[(q1 < n2) & (n2 <= q2)] <- i
          } else {
            nord2[n2 > qnorm(marg2[i])] <- i
          }
        }
        nord1 <- nord1 + 1
        nord2 <- nord2 + 1
        L_sigma[q, r] <- L_sigma[r, q] <-
          cor(nord1[order(nord1, decreasing = TRUE)], nord2[order(nord2)])
        U_sigma[q, r] <- U_sigma[r, q] <- cor(nord1[order(nord1)],
                                              nord2[order(nord2)])
      }
      if (q >= (k_cat + k_cont + 1) & q <= (k_cat + k_cont + k_pois) &
          r >= (k_cat + k_cont + k_pois + 1) &
          r <= (k_cat + k_cont + k_pois + k_nb)) {
        marg1 <- pois_marg[[q - (k_cat + k_cont)]]
        marg2 <- nb_marg[[r - (k_cat + k_cont + k_pois)]]
        n1 <- rnorms[, 1]
        n2 <- rnorms[, 2]
        nord1 <- numeric(length(n1))
        nord2 <- numeric(length(n2))
        for (i in 1:length(marg1)) {
          if (i != length(marg1)) {
            q1 <- qnorm(marg1[i])
            q2 <- qnorm(marg1[i + 1])
            nord1[(q1 < n1) & (n1 <= q2)] <- i
          } else {
            nord1[n1 > qnorm(marg1[i])] <- i
          }
        }
        for (i in 1:length(marg2)) {
          if (i != length(marg2)) {
            q1 <- qnorm(marg2[i])
            q2 <- qnorm(marg2[i + 1])
            nord2[(q1 < n2) & (n2 <= q2)] <- i
          } else {
            nord2[n2 > qnorm(marg2[i])] <- i
          }
        }
        nord1 <- nord1 + 1
        nord2 <- nord2 + 1
        L_sigma[q, r] <- L_sigma[r, q] <-
          cor(nord1[order(nord1, decreasing = TRUE)], nord2[order(nord2)])
        U_sigma[q, r] <- U_sigma[r, q] <- cor(nord1[order(nord1)],
                                              nord2[order(nord2)])
      }
      if (q >= (k_cat + k_cont + k_pois + 1) &
          q <= (k_cat + k_cont + k_pois + k_nb) &
          r >= (k_cat + k_cont + k_pois + 1) &
          r <= (k_cat + k_cont + k_pois + k_nb)) {
        marg1 <- nb_marg[[q - (k_cat + k_cont + k_pois)]]
        marg2 <- nb_marg[[r - (k_cat + k_cont + k_pois)]]
        n1 <- rnorms[, 1]
        n2 <- rnorms[, 2]
        nord1 <- numeric(length(n1))
        nord2 <- numeric(length(n2))
        for (i in 1:length(marg1)) {
          if (i != length(marg1)) {
            q1 <- qnorm(marg1[i])
            q2 <- qnorm(marg1[i + 1])
            nord1[(q1 < n1) & (n1 <= q2)] <- i
          } else {
            nord1[n1 > qnorm(marg1[i])] <- i
          }
        }
        for (i in 1:length(marg2)) {
          if (i != length(marg2)) {
            q1 <- qnorm(marg2[i])
            q2 <- qnorm(marg2[i + 1])
            nord2[(q1 < n2) & (n2 <= q2)] <- i
          } else {
            nord2[n2 > qnorm(marg2[i])] <- i
          }
        }
        nord1 <- nord1 + 1
        nord2 <- nord2 + 1
        L_sigma[q, r] <- L_sigma[r, q] <-
          cor(nord1[order(nord1, decreasing = TRUE)], nord2[order(nord2)])
        U_sigma[q, r] <- U_sigma[r, q] <- cor(nord1[order(nord1)],
                                              nord2[order(nord2)])
      }
    }
  }
  valid.state <- NULL
  if (!is.null(rho)) {
    valid.state <- TRUE
    for (i in 1:(k - 1)) {
      for (j in (i + 1):k) {
        if (rho[i, j] < L_sigma[i, j] | rho[i, j] > U_sigma[i, j]) {
          cat("Range error! Corr[", i, ",", j, "] must be between",
              round(L_sigma[i, j], 6), "and", round(U_sigma[i, j], 6), "\n")
          valid.state <- FALSE
        }
      }
    }
    if (valid.state == TRUE) cat("All correlations are in feasible range! \n")
    if (valid.state == FALSE)
      cat("Some correlations are not in feasible range! \n")
  }
  if (k_cont > 0) {
    return(list(rho = rho, L_rho = L_sigma, U_rho = U_sigma,
      constants = constants, sixth_correction = SixCorr,
      valid.pdf = Valid.PDF, valid.rho = valid.state))
  } else {
    return(list(rho = rho, L_rho = L_sigma, U_rho = U_sigma,
      valid.rho = valid.state))
  }
}
