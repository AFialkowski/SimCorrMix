#' @title Calculate Intermediate MVN Correlation for Ordinal - Poisson Variables: Correlation Method 1
#'
#' @description This function calculates a \code{k_cat x k_pois} intermediate matrix of correlations for the \code{k_cat} ordinal (\eqn{r >=
#'     2} categories) and \code{k_pois} Poisson variables required to produce the target correlations in \code{rho_cat_pois}. It extends the method of Amatya & Demirtas (2015, \doi{10.1080/00949655.2014.953534})
#'     to ordinal - Poisson pairs and allows for regular or zero-inflated Poisson variables.
#'     Here, the intermediate correlation between Z1 and Z2 (where Z1 is the standard normal variable discretized to produce an
#'     ordinal variable Y1, and Z2 is the standard normal variable used to generate a Poisson variable via the inverse CDF method) is
#'     calculated by dividing the target correlation by a correction factor.  The correction factor is the product of the
#'     upper Frechet-Hoeffding bound on the correlation between a Poisson variable and the normal variable used to generate it
#'     and a simulated GSC upper bound on the correlation between an ordinal variable and the normal variable used to generate it (see
#'     Demirtas & Hedeker, 2011, \doi{10.1198/tast.2011.10090}).  The function is used in \code{\link[SimCorrMix]{intercorr}} and
#'     \code{\link[SimCorrMix]{corrvar}}.  This function would not ordinarily be called by the user.
#'
#' @param rho_cat_pois a \code{k_cat x k_pois} matrix of target correlations among ordinal and Poisson variables; the Poisson variables
#'     should be ordered 1st regular, 2nd zero-inflated
#' @param marginal a list of length equal to \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param lam a vector of lambda (mean > 0) constants for the regular and zero-inflated Poisson variables (see \code{stats::dpois});
#'     the order should be 1st regular Poisson variables, 2nd zero-inflated Poisson variables
#' @param p_zip a vector of probabilities of structural zeros (not including zeros from the Poisson distribution) for the
#'     zero-inflated Poisson variables (see \code{VGAM::dzipois}); if \code{p_zip} = 0, \eqn{Y_{pois}} has a regular Poisson
#'     distribution; if \code{p_zip} is in (0, 1), \eqn{Y_{pois}} has a zero-inflated Poisson distribution;
#'     if \code{p_zip} is in \code{(-(exp(lam) - 1)^(-1), 0)}, \eqn{Y_{pois}} has a zero-deflated Poisson distribution and \code{p_zip}
#'     is not a probability; if \code{p_zip = -(exp(lam) - 1)^(-1)}, \eqn{Y_{pois}} has a positive-Poisson distribution
#'     (see \code{VGAM::dpospois}); if \code{length(p_zip) < length(lam)}, the missing values are set to 0 (and ordered 1st)
#' @param nrand the number of random numbers to generate in calculating the bound (default = 10000)
#' @param seed the seed used in random number generation (default = 1234)
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @importFrom VGAM qzipois
#' @export
#' @keywords correlation ordinal Poisson method1
#' @seealso \code{\link[SimCorrMix]{intercorr}}, \code{\link[SimCorrMix]{corrvar}}
#' @return a \code{k_cat x k_pois} matrix whose rows represent the \code{k_cat} ordinal variables and columns represent the \code{k_pois} Poisson variables
#' @references Amatya A & Demirtas H (2015). Simultaneous generation of multivariate mixed data with Poisson and normal marginals.
#'     Journal of Statistical Computation and Simulation, 85(15):3129-39. \doi{10.1080/00949655.2014.953534}.
#'
#' Demirtas H & Hedeker D (2011). A practical way for computing approximate lower and upper correlation bounds.
#'     American Statistician, 65(2):104-109. \doi{10.1198/tast.2011.10090}.
#'
#' Frechet M (1951). Sur les tableaux de correlation dont les marges sont donnees.  Ann. l'Univ. Lyon SectA, 14:53-77.
#'
#' Hoeffding W. Scale-invariant correlation theory. In: Fisher NI, Sen PK, editors. The collected works of Wassily Hoeffding.
#'     New York: Springer-Verlag; 1994. p. 57-107.
#'
#' Yahav I & Shmueli G (2012). On Generating Multivariate Poisson Data in Management Science Applications. Applied Stochastic
#'     Models in Business and Industry, 28(1):91-102. \doi{10.1002/asmb.901}.
#'
#' Yee TW (2018). VGAM: Vector Generalized Linear and Additive Models. R package version 1.0-5. \url{https://CRAN.R-project.org/package=VGAM}.
#'
intercorr_cat_pois <- function(rho_cat_pois = NULL, marginal = list(),
                               lam = NULL, p_zip = 0, nrand = 100000,
                               seed = 1234) {
  Sigma_cat_pois <- matrix(1, nrow = nrow(rho_cat_pois),
                           ncol = ncol(rho_cat_pois))
  if (length(p_zip) < length(lam))
    p_zip <- c(rep(0, length(lam) - length(p_zip)), p_zip)
  set.seed(seed)
  u <- runif(nrand, 0, 1)
  set.seed(seed)
  n2 <- rnorm(nrand, 0, 1)
  for (i in 1:nrow(rho_cat_pois)) {
    for (j in 1:ncol(rho_cat_pois)) {
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
      chat_pois <- cor(qzipois(u, lam[j], pstr0 = p_zip[j]), qnorm(u, 0, 1))
      Sigma_cat_pois[i, j] <-
        rho_cat_pois[i, j]/(chat_pois * cor(yord[order(yord)], n2[order(n2)]))
      if (Sigma_cat_pois[i, j] > 1) Sigma_cat_pois[i, j] <- 1
      if (Sigma_cat_pois[i, j] < -1) Sigma_cat_pois[i, j] <- -1
    }
  }
  return(Sigma_cat_pois)
}
