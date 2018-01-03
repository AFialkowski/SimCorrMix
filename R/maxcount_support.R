#' @title Calculate Maximum Support Value for Count Variables: Correlation Method 2
#'
#' @description This function calculates the maximum support value for count variables by extending the method of Barbiero &
#'     Ferrari (2015, \doi{10.1002/asmb.2072}) to include regular and zero-inflated Poisson and Negative Binomial variables.  In order for
#'     count variables to be treated as ordinal in the calculation of the intermediate MVN correlation matrix, their infinite support must
#'     be truncated (made finite).  This is done by setting the total cumulative probability equal to 1 - a small user-specified value
#'     (\code{pois_eps} or \code{nb_eps}).  The maximum support value equals the inverse CDF applied to this result.  The truncation values
#'     may differ for each variable.  The function is used in \code{\link[SimCorrMix]{intercorr2}} and \code{\link[SimCorrMix]{corrvar2}} and
#'     would not ordinarily be called by the user.
#'
#' @param k_pois the number of Poisson variables
#' @param k_nb the number of Negative Binomial variables
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
#' @param prob a vector of success probability parameters for the NB variables; order the same as in \code{size}
#' @param mu a vector of mean parameters for the NB variables (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL); order the same as in \code{size}; for zero-inflated NB this refers to
#'     the mean of the NB distribution (see \code{\link[VGAM;Zinegbin]{dzinegbin}})
#' @param p_zinb a vector of probabilities of structural zeros (not including zeros from the NB distribution) for the zero-inflated NB variables
#'     (see \code{\link[VGAM;Zinegbin]{dzinegbin}}); if \code{p_zinb} = 0, \eqn{Y_{nb}} has a regular NB distribution;
#'     if \code{p_zinb} is in \code{(-prob^size/(1 - prob^size),} \code{0)}, \eqn{Y_{nb}} has a zero-deflated NB distribution and \code{p_zinb}
#'     is not a probability; if \code{p_zinb = -prob^size/(1 - prob^size)}, \eqn{Y_{nb}} has a positive-NB distribution (see
#'     \code{\link[VGAM;Posnegbin]{dposnegbin}}); if \code{length(p_zinb) < length(size)}, the missing values are set to 0 (and ordered 1st)
#' @param pois_eps a vector of length \code{k_pois} containing total cumulative probability truncation values; if none are provided,
#'     the default is 0.0001 for each variable
#' @param nb_eps a vector of length \code{k_nb} containing total cumulative probability truncation values; if none are provided,
#'     the default is 0.0001 for each variable
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @importFrom VGAM qzipois qzinegbin
#' @export
#' @keywords Poisson NegativeBinomial method2
#' @seealso \code{\link[SimCorrMix]{intercorr2}}, \code{\link[SimCorrMix]{corrvar2}}
#' @return a data.frame with \code{k_pois + k_nb} rows; the column names are:
#' @return \code{Distribution} Poisson or Negative Binomial
#' @return \code{Number} the variable index
#' @return \code{Max} the maximum support value
#' @references
#' Barbiero A & Ferrari PA (2015). Simulation of correlated Poisson variables. Applied Stochastic Models in
#'     Business and Industry, 31:669-80. \doi{10.1002/asmb.2072}.
#'
#' Ferrari PA, Barbiero A (2012). Simulating ordinal data, Multivariate Behavioral Research, 47(4):566-589. \doi{10.1080/00273171.2012.692630}.
#'
#' Yee TW (2017). VGAM: Vector Generalized Linear and Additive Models. \cr \url{https://CRAN.R-project.org/package=VGAM}.
#'
maxcount_support <- function(k_pois = 0, k_nb = 0, lam = NULL, p_zip = 0,
                             size = NULL, prob = NULL, mu = NULL, p_zinb = 0,
                             pois_eps = NULL, nb_eps = NULL) {
  max_support <- matrix(1, nrow = k_pois + k_nb, ncol = 2)
  if (length(p_zinb) < length(size))
    p_zinb <- c(rep(0, length(size) - length(p_zinb)), p_zinb)
  if (length(p_zip) < length(lam))
    p_zip <- c(rep(0, length(lam) - length(p_zip)), p_zip)
  if (k_pois > 0) {
    for (i in 1:k_pois) {
      max_support[i, ] <- append(i, qzipois(1 - pois_eps[i], lam[i],
                                            pstr0 = p_zip[i]))
    }
  }
  if (k_nb > 0) {
    if (length(prob) > 0)
      mu <- size * (1 - prob)/prob
    for (i in (k_pois + 1):(k_pois + k_nb)) {
      max_support[i, ] <- append(i, qzinegbin(1 - nb_eps[i - k_pois],
        size = size[i - k_pois], munb = mu[i - k_pois],
        pstr0 = p_zinb[i - k_pois]))
    }
  }
  max_support <- cbind(append(rep("Poisson", k_pois),
                              rep("Neg_Bin", k_nb)),
                       as.data.frame(max_support))
  colnames(max_support) <- c("Distribution", "Number", "Max")
  return(max_support)
}
