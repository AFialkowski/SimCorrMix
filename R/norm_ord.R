#' @title Calculate Correlations of Ordinal Variables Obtained from Discretizing Normal Variables
#'
#' @description This function calculates the correlation of ordinal variables (or variables treated as "ordinal"), with given marginal
#'     distributions, obtained from discretizing standard normal variables with a specified correlation matrix.  The function modifies
#'     Barbiero & Ferrari's \code{\link[GenOrd]{contord}} function in \code{\link[GenOrd]{GenOrd-package}}.  It uses
#'     \code{\link[mvtnorm]{pmvnorm}} function from the \strong{mvtnorm} package to calculate multivariate normal cumulative probabilities
#'     defined by the normal quantiles obtained at \code{marginal} and the supplied correlation matrix \code{Sigma}.  This function is used
#'     within \code{\link[SimCorrMix]{ord_norm}} and would not ordinarily be called by the user.
#'
#' @param marginal a list of length equal to the number of variables; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1)
#' @param Sigma the correlation matrix of the multivariate standard normal variable
#' @param support a list of length equal to the number of variables; the i-th element is a vector of containing the r
#'     ordered support values; if not provided (i.e. support = list()), the default is for the i-th element to be the vector 1, ..., r
#' @param Spearman if TRUE, Spearman's correlations are used (and support is not required); if FALSE (default) Pearson's correlations
#'     are used
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @export
#' @keywords correlation ordinal continuous
#' @seealso \code{\link[SimCorrMix]{ord_norm}}
#' @return the correlation matrix of the ordinal variables
#' @references
#' Please see references in \code{\link[SimCorrMix]{ord_norm}}.
#'
#' Genz A, Bretz F, Miwa T, Mi X, Leisch F, Scheipl F, Hothorn T (2017).  mvtnorm: Multivariate Normal and t Distributions.
#'     R package version 1.0-6.  \url{https://CRAN.R-project.org/package=mvtnorm}
#'
#' Genz A, Bretz F (2009), Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195.,
#'     Springer-Verlag, Heidelberg. ISBN 978-3-642-01688-2
#'
norm_ord <- function(marginal = list(), Sigma = NULL, support = list(),
                     Spearman = FALSE) {
  if (!all(unlist(lapply(marginal,
                         function(x) (sort(x) == x & min(x) > 0 &
                                      max(x) < 1))))) {
    stop("Error in given marginal distributions!")
  }
  if (!isSymmetric(Sigma) |
      !all(diag(Sigma) == 1)) {
    stop("Correlation matrix not valid!")
  }
  k_cat <- length(marginal)
  Sigmaord <- diag(1, k_cat)
  if (length(support) == 0) {
    for (i in 1:k_cat) {
      if (Spearman == TRUE) {
        s1 <- c(marginal[[i]], 1)
        s2 <- c(0, marginal[[i]])
        support[[i]] <- (s1 + s2)/2
      } else
        support[[i]] <- 1:(length(marginal[[i]]) + 1)
    }
  }
  L <- vector("list", k_cat)
  k.cat <- numeric(k_cat)
  for (i in 1:k_cat) {
    k.cat[i] <- length(marginal[[i]]) + 1
    L[[i]] <- qnorm(marginal[[i]])
    L[[i]] <- c(-Inf, L[[i]], Inf)
  }
  for (q in 1:(k_cat - 1)) {
    for (r in (q + 1):k_cat) {
      p.qr <- matrix(0, k.cat[q], k.cat[r])
      for (i in 1:k.cat[q]) {
        for (j in 1:k.cat[r]) {
          low <- rep(-Inf, k_cat)
          upp <- rep(Inf, k_cat)
          low[q] <- L[[q]][i]
          low[r] <- L[[r]][j]
          upp[q] <- L[[q]][i + 1]
          upp[r] <- L[[r]][j + 1]
          p.qr[i, j] <- pmvnorm(low, upp, rep(0, k_cat), corr = Sigma)
        }
      }
      mean.y <- sum(apply(p.qr, 2, sum) * support[[r]])
      sigma.y <- sqrt(sum(apply(p.qr, 2, sum) * support[[r]]^2) - mean.y^2)
      mean.x <- sum(apply(p.qr, 1, sum) * support[[q]])
      sigma.x <- sqrt(sum(apply(p.qr, 1, sum) * support[[q]]^2) - mean.x^2)
      mean.qr <- support[[q]] %*% t(support[[r]])
      mu.qr <- sum(mean.qr * p.qr)
      cov.xy <- mu.qr - mean.x * mean.y
      cor.xy <- cov.xy/(sigma.x * sigma.y)
      Sigmaord[q, r] <- cor.xy
      Sigmaord[r, q] <- Sigmaord[q, r]
    }
  }
  Sigmaord
}
