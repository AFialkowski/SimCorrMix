#' @title Calculate Intermediate MVN Correlation for Ordinal - Negative Binomial Variables: Correlation Method 1
#'
#' @description This function calculates the \code{k_cat x k_nb} intermediate matrix of correlations for the \code{k_cat} ordinal (\eqn{r >=
#'     2} categories) and \code{k_nb} Negative Binomial variables required to produce the target correlations in \code{rho_cat_nb}. It extends the method of Amatya & Demirtas (2015, \doi{10.1080/00949655.2014.953534})
#'     to ordinal - Negative Binomial pairs and allows for regular or zero-inflated NB variables.  Here, the intermediate correlation between Z1 and Z2 (where Z1 is the standard normal variable
#'     discretized to produce an ordinal variable Y1, and Z2 is the standard normal variable used to generate a Negative Binomial
#'     variable via the inverse CDF method) is calculated by dividing the target correlation by a correction factor.  The
#'     correction factor is the product of the upper Frechet-Hoeffding bound on the correlation between a Negative Binomial variable
#'     and the normal variable used to generate it and a simulated GSC upper bound on the correlation between an ordinal variable and the normal variable used to generate it (see Demirtas & Hedeker, 2011,
#'     \doi{10.1198/tast.2011.10090}).  The function is used in \code{\link[SimCorrMix]{intercorr}} and \code{\link[SimCorrMix]{corrvar}}.
#'     This function would not ordinarily be called by the user.
#'
#' @param rho_cat_nb a \code{k_cat x k_nb} matrix of target correlations among ordinal and Negative Binomial variables; the NB variables
#'     should be ordered 1st regular, 2nd zero-inflated
#' @param marginal a list of length equal to \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{stats::dnbinom}); the order should be
#'     1st regular NB variables, 2nd zero-inflated NB variables
#' @param mu a vector of mean parameters for the NB variables (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL); order the same as in \code{size}; for zero-inflated NB this refers to
#'     the mean of the NB distribution (see \code{VGAM::dzinegbin})
#' @param p_zinb a vector of probabilities of structural zeros (not including zeros from the NB distribution) for the zero-inflated NB variables
#'     (see \code{VGAM::dzinegbin}); if \code{p_zinb} = 0, \eqn{Y_{nb}} has a regular NB distribution;
#'     if \code{p_zinb} is in \code{(-prob^size/(1 - prob^size),} \code{0)}, \eqn{Y_{nb}} has a zero-deflated NB distribution and \code{p_zinb}
#'     is not a probability; if \code{p_zinb = -prob^size/(1 - prob^size)}, \eqn{Y_{nb}} has a positive-NB distribution (see
#'     \code{VGAM::dposnegbin}); if \code{length(p_zinb) < length(size)}, the missing values are set to 0 (and ordered 1st)
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @importFrom VGAM qzinegbin
#' @export
#' @keywords correlation ordinal NegativeBinomial method1
#' @seealso \code{\link[SimCorrMix]{intercorr}}, \code{\link[SimCorrMix]{corrvar}}
#' @return a \code{k_cat x k_nb} matrix whose rows represent the \code{k_cat} ordinal variables and columns represent the
#'     \code{k_nb} Negative Binomial variables
#' @references Please see references for \code{\link[SimCorrMix]{intercorr_cat_pois}}.
#'
intercorr_cat_nb <- function(rho_cat_nb = NULL, marginal = list(),
                             size = NULL, mu = NULL, p_zinb = 0,
                             nrand = 100000, seed = 1234) {
  Sigma_cat_nb <- matrix(1, nrow = nrow(rho_cat_nb), ncol = ncol(rho_cat_nb))
  if (length(p_zinb) < length(size))
    p_zinb <- c(rep(0, length(size) - length(p_zinb)), p_zinb)
  set.seed(seed)
  u <- runif(nrand, 0, 1)
  set.seed(seed)
  n2 <- rnorm(nrand, 0, 1)
  for (i in 1:nrow(rho_cat_nb)) {
    for (j in 1:ncol(rho_cat_nb)) {
      yord <- numeric(length(n2))
      for (r in 1:length(marginal[[i]])) {
        if (r != length(marginal[[i]])) {
          q1 <- qnorm(marginal[[i]][r])
          q2 <- qnorm(marginal[[i]][r + 1])
          yord[(q1 < n2) & (n2 <= q2)] <- r
        } else {
          yord[n2 > qnorm(marginal[[i]][r])] <- r
        }
      }
      yord <- yord + 1
      chat_nb <- cor(qzinegbin(u, size = size[j], munb = mu[j], pstr0 = p_zinb[j]),
                     qnorm(u, 0, 1))
      Sigma_cat_nb[i, j] <-
        rho_cat_nb[i, j]/(chat_nb * cor(yord[order(yord)], n2[order(n2)]))
      if (Sigma_cat_nb[i, j] > 1) Sigma_cat_nb[i, j] <- 1
      if (Sigma_cat_nb[i, j] < -1) Sigma_cat_nb[i, j] <- -1
    }
  }
  return(Sigma_cat_nb)
}
