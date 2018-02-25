#' @title Approximate Correlation between Two Continuous Mixture Variables M1 and M2
#'
#' @description This function approximates the expected correlation between two continuous mixture variables \eqn{M1} and \eqn{M2} based on
#'     their mixing proportions, component means, component standard deviations, and correlations between components across variables.
#'     The equations can be found in the \strong{Expected Cumulants and Correlations for Continuous Mixture Variables} vignette.  This
#'     function can be used to see what combination of component correlations gives a desired correlation between \eqn{M1} and \eqn{M2}.
#'
#' @param mix_pis a list of length 2 with 1st component a vector of mixing probabilities that sum to 1 for component distributions of
#'     \eqn{M1} and likewise for 2nd component and \eqn{M2}
#' @param mix_mus a list of length 2 with 1st component a vector of means for component distributions of \eqn{M1} and likewise for 2nd
#'     component and \eqn{M2}
#' @param mix_sigmas a list of length 2 with 1st component a vector of standard deviations for component distributions of \eqn{M1} and
#'     likewise for 2nd component and \eqn{M2}
#' @param p_M1M2 a matrix of correlations with rows corresponding to \eqn{M1} and columns corresponding to \eqn{M2}; i.e.,
#'     \code{p_M1M2[1, 2]} is the correlation between the 1st component of \eqn{M1} and the 2nd component of \eqn{M2}
#'
#' @export
#' @seealso \code{\link[SimCorrMix]{rho_M1Y}}
#' @keywords mixture correlation
#' @return the expected correlation between M1 and M2
#' @references
#' Davenport JW, Bezder JC, & Hathaway RJ (1988). Parameter Estimation for Finite Mixture Distributions.
#'     Computers & Mathematics with Applications, 15(10):819-28.
#'
#' Pearson RK (2011). Exploring Data in Engineering, the Sciences, and Medicine. In. New York: Oxford University Press.
#'
#' @examples
#' # M1 is mixture of N(-2, 1) and N(2, 1);
#' # M2 is mixture of Logistic(0, 1), Chisq(4), and Beta(4, 1.5)
#' # pairwise correlation between components across M1 and M2 set to 0.35
#' L <- calc_theory("Logistic", c(0, 1))
#' C <- calc_theory("Chisq", 4)
#' B <- calc_theory("Beta", c(4, 1.5))
#' mix_pis <- list(c(0.4, 0.6), c(0.3, 0.2, 0.5))
#' mix_mus <- list(c(-2, 2), c(L[1], C[1], B[1]))
#' mix_sigmas <- list(c(1, 1), c(L[2], C[2], B[2]))
#' p_M11M21 <- p_M11M22 <- p_M11M23 <- 0.35
#' p_M12M21 <- p_M12M22 <- p_M12M23 <- 0.35
#' p_M1M2 <- matrix(c(p_M11M21, p_M11M22, p_M11M23, p_M12M21, p_M12M22,
#'   p_M12M23), 2, 3, byrow = TRUE)
#' rhoM1M2 <- rho_M1M2(mix_pis, mix_mus, mix_sigmas, p_M1M2)
#' rhoM1M2
#'
rho_M1M2 <- function(mix_pis = list(), mix_mus = list(), mix_sigmas = list(),
                     p_M1M2 = NULL) {
  if (all(dim(p_M1M2) %in% length(mix_pis[[1]]))) {
    if (all(p_M1M2 == diag(length(mix_pis[[1]])))) return(1)
  }
  one <- mix_pis[[1]] * mix_sigmas[[1]]
  two <- apply(p_M1M2, 1, function(x) sum(mix_pis[[2]] * mix_sigmas[[2]] * x))
  Var <- mapply(function(x, y, z) sum(x * (z^2 + y^2)) - (sum(x * y))^2,
                mix_pis, mix_mus, mix_sigmas)
  rhoM1M2 <- sum(one * two)/sqrt(Var[1] * Var[2])
  rhoM1M2
}
