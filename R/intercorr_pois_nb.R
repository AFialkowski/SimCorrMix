#' @title Calculate Intermediate MVN Correlation for Poisson - Negative Binomial Variables: Correlation Method 1
#'
#' @description This function calculates a \code{k_pois x k_nb} intermediate matrix of correlations for the
#'     Poisson and Negative Binomial variables by extending the method of Yahav & Shmueli (2012, \doi{10.1002/asmb.901}). The intermediate
#'     correlation between Z1 and Z2 (the standard normal variables used to generate the Poisson and Negative Binomial variables Y1 and Y2
#'     via the inverse CDF method) is calculated using a logarithmic transformation of the target correlation.  First, the upper and lower
#'     Frechet-Hoeffding bounds (mincor, maxcor) on \eqn{\rho_{Y1, Y2}} are simulated.  Then the intermediate correlation is found as follows:
#'     \deqn{\rho_{Z1, Z2} = \frac{1}{b} * log(\frac{\rho_{Y1, Y2} - c}{a}),}
#'     where \eqn{a = -(maxcor * mincor)/(maxcor + mincor)}, \eqn{b = log((maxcor + a)/a)}, and \eqn{c = -a}.
#'     The function adapts code from Amatya & Demirtas' (2016) package \code{\link[PoisNor]{PoisNor-package}} by:
#'
#'     1) allowing specifications for the number of random variates and the seed for reproducibility
#'
#'     2) providing the following checks: if \code{Sigma_(Z1, Z2)} > 1, \code{Sigma_(Z1, Z2)} is set to 1; if \code{Sigma_(Z1, Z2)} < -1,
#'     \code{Sigma_(Z1, Z2)} is set to -1
#'
#'     3) simulating regular and zero-inflated Poisson and Negative Binomial variables.
#'
#'     The function is used in \code{\link[SimCorrMix]{intercorr}} and \code{\link[SimCorrMix]{corrvar}} and would not ordinarily be called by the user.
#'
#' @param rho_pois_nb a \code{k_pois x k_nb} matrix of target correlations; order of each type should be 1st regular, 2nd zero-inflated
#' @param lam a vector of lambda (mean > 0) constants for the regular and zero-inflated Poisson variables (see \code{\link[stats;Poisson]{dpois}});
#'     the order should be 1st regular Poisson variables, 2nd zero-inflated Poisson variables
#' @param p_zip a vector of probabilities of structural zeros (not including zeros from the Poisson distribution) for the
#'     zero-inflated Poisson variables (see \code{\link[VGAM;Zipois]{dzipois}}); if \code{p_zip} = 0, \eqn{Y_{pois}} has a regular Poisson
#'     distribution; if \code{p_zip} is in (0, 1), \eqn{Y_{pois}} has a zero-inflated Poisson distribution;
#'     if \code{p_zip} is in \code{(-(exp(lam) - 1)^(-1), 0)}, \eqn{Y_{pois}} has a zero-deflated Poisson distribution and \code{p_zip}
#'     is not a probability; if \code{p_zip = -(exp(lam) - 1)^(-1)}, \eqn{Y_{pois}} has a positive-Poisson distribution
#'     (see \code{\link[VGAM;Pospois]{dpospois}}); if \code{length(p_zip) < length(lam)}, the missing values are set to 0 (and ordered 1st)
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats;NegBinomial]{dnbinom}}); the order should be
#'     1st regular NB variables, 2nd zero-inflated NB variables
#' @param mu a vector of mean parameters for the NB variables (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL); order the same as in \code{size}; for zero-inflated NB this refers to
#'     the mean of the NB distribution (see \code{\link[VGAM;Zinegbin]{dzinegbin}})
#' @param p_zinb a vector of probabilities of structural zeros (not including zeros from the NB distribution) for the zero-inflated NB variables
#'     (see \code{\link[VGAM;Zinegbin]{dzinegbin}}); if \code{p_zinb} = 0, \eqn{Y_{nb}} has a regular NB distribution;
#'     if \code{p_zinb} is in \code{(-prob^size/(1 - prob^size),} \code{0)}, \eqn{Y_{nb}} has a zero-deflated NB distribution and \code{p_zinb}
#'     is not a probability; if \code{p_zinb = -prob^size/(1 - prob^size)}, \eqn{Y_{nb}} has a positive-NB distribution (see
#'     \code{\link[VGAM;Posnegbin]{dposnegbin}}); if \code{length(p_zinb) < length(size)}, the missing values are set to 0 (and ordered 1st)
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @importFrom VGAM qzipois qzinegbin
#' @export
#' @keywords correlation Poisson NegativeBinomial method1
#' @seealso \code{\link[SimCorrMix]{intercorr_pois}}, \code{\link[SimCorrMix]{intercorr_nb}},
#'     \code{\link[SimCorrMix]{intercorr}}, \code{\link[SimCorrMix]{corrvar}}
#' @return the \code{k_pois x k_nb} intermediate correlation matrix whose rows represent the \code{k_pois} Poisson variables and
#'     columns represent the \code{k_nb} Negative Binomial variables
#' @references Please see references for \code{\link[SimCorrMix]{intercorr_pois}}.
#'
intercorr_pois_nb <- function(rho_pois_nb = NULL, lam = NULL, p_zip = 0,
                              size = NULL, mu = NULL, p_zinb = 0,
                              nrand = 100000, seed = 1234) {
  if (length(p_zinb) < length(size))
    p_zinb <- c(rep(0, length(size) - length(p_zinb)), p_zinb)
  if (length(p_zip) < length(lam))
    p_zip <- c(rep(0, length(lam) - length(p_zip)), p_zip)
  Sigma_pois_nb <- matrix(1, nrow(rho_pois_nb), ncol(rho_pois_nb))
  set.seed(seed)
  u <- runif(nrand, 0, 1)
  for (i in 1:nrow(rho_pois_nb)) {
    for (j in 1:ncol(rho_pois_nb)) {
      maxcor <- cor(qzipois(u, lam[i], pstr0 = p_zip[i]),
                    qzinegbin(u, size = size[j], munb = mu[j],
                              pstr0 = p_zinb[j]))
      mincor <- cor(qzipois(u, lam[i], pstr0 = p_zip[i]),
                    qzinegbin(1 - u, size = size[j], munb = mu[j],
                              pstr0 = p_zinb[j]))
      a <- -(maxcor * mincor)/(maxcor + mincor)
      b <- log((maxcor + a)/a)
      c <- -a
      Sigma_pois_nb[i, j] <- (1/b) * log((rho_pois_nb[i, j] - c)/a)
      if (Sigma_pois_nb[i, j] > 1) Sigma_pois_nb[i, j] <- 1
      if (Sigma_pois_nb[i, j] < -1) Sigma_pois_nb[i, j] <- 1
    }
  }
  return(Sigma_pois_nb)
}
