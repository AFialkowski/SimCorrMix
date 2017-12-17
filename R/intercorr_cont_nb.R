#' @title Calculate Intermediate MVN Correlation for Continuous - Negative Binomial Variables: Correlation Method 1
#'
#' @description This function calculates a \code{k_cont x k_nb} intermediate matrix of correlations for the \code{k_cont} continuous and
#'     \code{k_nb} Negative Binomial variables. It extends the method of Amatya & Demirtas (2015, \doi{10.1080/00949655.2014.953534}) to
#'     continuous variables generated using Headrick's fifth-order polynomial transformation and regular or zero-inflated NB variables.
#'     Here, the intermediate correlation between Z1 and Z2 (where Z1 is the standard normal variable transformed using Headrick's fifth-order
#'     or Fleishman's third-order method to produce a continuous variable Y1, and Z2 is the standard normal variable used to generate a
#'     Negative Binomial variable via the inverse CDF method) is calculated by dividing the target correlation by a correction factor.
#'     The correction factor is the product of the upper Frechet-Hoeffding bound on the correlation between a Negative Binomial variable and
#'     the normal variable used to generate it and the power method correlation (described in Headrick & Kowalchuk, 2007,
#'     \doi{10.1080/10629360600605065}) between Y1 and Z1.  The function is used in \code{\link[SimCorrMix]{intercorr}} and
#'     \code{\link[SimCorrMix]{corrvar}}.  This function would not ordinarily be called by the user.
#'
#' @param method the method used to generate the \code{k_cont} continuous variables.  "Fleishman" uses a third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param constants a matrix with \code{k_cont} rows, each a vector of constants c0, c1, c2, c3 (if \code{method} = "Fleishman") or
#'     c0, c1, c2, c3, c4, c5 (if \code{method} = "Polynomial"), like that returned by \code{\link[SimMultiCorrData]{find_constants}}
#' @param rho_cont_nb a \code{k_cont x k_nb} matrix of target correlations among continuous and Negative Binomial variables; the NB variables
#'     should be ordered 1st regular, 2nd zero-inflated
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{dnbinom}}); the order should be
#'     1st regular NB variables, 2nd zero-inflated NB variables
#' @param mu a vector of mean parameters for the NB variables; order the same as in \code{size}; for zero-inflated NB this refers to
#'     the mean of the NB distribution (see \code{\link[VGAM]{dzinegbin}})
#' @param p_zinb a vector of probabilities of structural zeros (not including zeros from the NB distribution) for the zero-inflated NB variables
#'     (see \code{\link[VGAM]{dzinegbin}}); if \code{p_zinb} = 0, \eqn{Y_{nb}} has a regular NB distribution;
#'     if \code{p_zinb} is in \code{(-prob^size/(1 - prob^size),} \code{0)}, \eqn{Y_{nb}} has a zero-deflated NB distribution and \code{p_zinb}
#'     is not a probability; if \code{p_zinb = -prob^size/(1 - prob^size)}, \eqn{Y_{nb}} has a positive-NB distribution (see
#'     \code{\link[VGAM]{dposnegbin}}); if \code{length(p_zinb) < length(size)}, the missing values are set to 0 (and ordered 1st)
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @importFrom VGAM qzinegbin
#' @export
#' @keywords correlation continuous NegativeBinomial method1
#' @seealso \code{\link[SimMultiCorrData]{find_constants}},
#'     \code{\link[SimCorrMix]{intercorr}}, \code{\link[SimCorrMix]{corrvar}}
#' @return a \code{k_cont x k_nb} matrix whose rows represent the \code{k_cont} continuous variables and columns represent the
#'     \code{k_nb} Negative Binomial variables
#' @references Please see references for \code{\link[SimCorrMix]{intercorr_cont_pois}}.
#'
intercorr_cont_nb <- function(method = c("Fleishman", "Polynomial"),
                              constants = NULL, rho_cont_nb = NULL,
                              size = NULL, mu = NULL, p_zinb = 0,
                              nrand = 100000, seed = 1234) {
  Sigma_cont_nb <- matrix(1, nrow = nrow(rho_cont_nb),
                          ncol = ncol(rho_cont_nb))
  if (length(p_zinb) < length(size))
    p_zinb <- c(rep(0, length(size) - length(p_zinb)), p_zinb)
  set.seed(seed)
  u <- runif(nrand, 0, 1)
  for (i in 1:nrow(rho_cont_nb)) {
    for (j in 1:ncol(rho_cont_nb)) {
      chat_nb <- cor(qzinegbin(u, size = size[j], munb = mu[j],
        pstr0 = p_zinb[j]), qnorm(u, 0, 1))
      Sigma_cont_nb[i, j] <-
        rho_cont_nb[i, j]/(chat_nb * power_norm_corr(constants[i, ], method))
      if (Sigma_cont_nb[i, j] > 1) Sigma_cont_nb[i, j] <- 1
      if (Sigma_cont_nb[i, j] < -1) Sigma_cont_nb[i, j] <- -1
    }
  }
  return(Sigma_cont_nb)
}
