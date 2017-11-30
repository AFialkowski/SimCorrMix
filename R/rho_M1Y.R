#' @title Theoretical Correlation between Continuous Mixture Variable M1 and Random Variable Y
#'
#' @description This function calculates the expected correlation between a continous mixture variables \eqn{M1} and another random
#'     variable \eqn{Y} based on the mixing proportions, component means, and component standard deviations of \eqn{M1} and correlations
#'     between components of \eqn{M1} and \eqn{Y}.  The equations can be found in the \strong{Expected Cumulants and Correlations for
#'     Continuous Mixture Variables} vignette.  This function can be used to see what combination of correlations between components
#'     of \eqn{M1} and \eqn{Y} gives a desired correlation between \eqn{M1} and \eqn{Y}.
#'
#' @param mix_pis a vector of mixing probabilities that sum to 1 for component distributions of \eqn{M1}
#' @param mix_mus a vector of means for component distributions of \eqn{M1}
#' @param mix_sigmas a vector of standard deviations for component distributions of \eqn{M1}
#' @param p_M1Y a vector of correlations between the components of \eqn{M1} and \eqn{Y}; i.e.,
#'     \code{p_M1Y[1]} is the correlation between the 1st component of \eqn{M1} and \eqn{Y}
#'
#' @export
#' @seealso \code{\link[SimCorrMix]{rho_M1Y}}
#' @keywords mixture correlation
#' @references
#' Please see references for \code{\link[SimCorrMix]{rho_M1Y}}.
#'
#' @examples
#' # M1 is mixture of N(-2, 1) and N(2, 1); C1 is a continuous non-mixture
#' # variable (but could also be an ordinal or count variable)
#' # pairwise correlation between components of M1 and C1 set to 0.35
#' L <- calc_theory("Logistic", c(0, 1))
#' C <- calc_theory("Chisq", 4)
#' B <- calc_theory("Beta", c(4, 1.5))
#' mix_pis <- list(c(0.4, 0.6), c(0.3, 0.2, 0.5))
#' mix_mus <- list(c(-2, 2), c(L[1], C[1], B[1]))
#' mix_sigmas <- list(c(1, 1), c(L[2], C[2], B[2]))
#' p_M11C1 <- p_M12C1 <- 0.35
#' p_M1C1 <- c(p_M11C1, p_M12C1)
#' rho_M1C1 <- rho_M1Y(mix_pis[[1]], mix_mus[[1]], mix_sigmas[[1]], p_M1C1)
#' rho_M1C1
#'
rho_M1Y <- function(mix_pis = NULL, mix_mus = NULL, mix_sigmas = NULL,
                    p_M1Y = NULL) {
  Var <- sum(mix_pis * (mix_sigmas^2 + mix_mus^2)) - (sum(mix_pis * mix_mus))^2
  rhoM1Y <- sum(mix_pis * mix_sigmas * p_M1Y)/sqrt(Var)
  rhoM1Y
}
