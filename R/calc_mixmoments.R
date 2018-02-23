#' @title Find Standardized Cumulants of a Continuous Mixture Distribution by Method of Moments
#'
#' @description This function uses the method of moments to calculate the expected mean, standard deviation, skewness,
#'     standardized kurtosis, and standardized fifth and sixth cumulants for a continuous mixture variable based on the distributions
#'     of its components.  The result can be used as input to \code{\link[SimMultiCorrData]{find_constants}} or for comparison to a
#'     simulated mixture variable from \code{\link[SimCorrMix]{contmixvar1}}, \code{\link[SimCorrMix]{corrvar}}, or
#'     \code{\link[SimCorrMix]{corrvar2}}.  See the \bold{Expected Cumulants and Correlations for Continuous Mixture Variables} vignette
#'     for equations of the cumulants.
#'
#' @param mix_pis a vector of mixing probabilities that sum to 1 for the component distributions
#' @param mix_mus a vector of means for the component distributions
#' @param mix_sigmas a vector of standard deviations for the component distributions
#' @param mix_skews a vector of skew values for the component distributions
#' @param mix_skurts a vector of standardized kurtoses for the component distributions
#' @param mix_fifths a vector of standardized fifth cumulants for the component distributions; keep NULL if using \code{method} = "Fleishman"
#'     to generate continuous variables
#' @param mix_sixths a vector of standardized sixth cumulants for the component distributions; keep NULL if using \code{method} = "Fleishman"
#'     to generate continuous variables
#'
#' @export
#' @keywords cumulants mixture
#' @return A vector of the mean, standard deviation, skewness, standardized kurtosis, and standardized fifth and sixth cumulants
#' @references Please see references for \code{\link[SimCorrMix]{SimCorrMix}}.
#'
#' @examples
#' # Two mixture variables: 1st is mixture of Normal(-2, 1) and Normal(2, 1),
#' # 2nd is mixture of Logistic(0, 1), Chisq(4), and Beta(4, 1.5)
#' L <- calc_theory("Logistic", c(0, 1))
#' C <- calc_theory("Chisq", 4)
#' B <- calc_theory("Beta", c(4, 1.5))
#' mix_pis <- list(c(0.4, 0.6), c(0.3, 0.2, 0.5))
#' mix_mus <- list(c(-2, 2), c(L[1], C[1], B[1]))
#' mix_sigmas <- list(c(1, 1), c(L[2], C[2], B[2]))
#' mix_skews <- list(rep(0, 2), c(L[3], C[3], B[3]))
#' mix_skurts <- list(rep(0, 2), c(L[4], C[4], B[4]))
#' mix_fifths <- list(rep(0, 2), c(L[5], C[5], B[5]))
#' mix_sixths <- list(rep(0, 2), c(L[6], C[6], B[6]))
#' Nstcum <- calc_mixmoments(mix_pis[[1]], mix_mus[[1]], mix_sigmas[[1]],
#'   mix_skews[[1]], mix_skurts[[1]], mix_fifths[[1]], mix_sixths[[1]])
#' Nstcum
#' Mstcum <- calc_mixmoments(mix_pis[[2]], mix_mus[[2]], mix_sigmas[[2]],
#'   mix_skews[[2]], mix_skurts[[2]], mix_fifths[[2]], mix_sixths[[2]])
#' Mstcum
#'
calc_mixmoments <- function(mix_pis = NULL, mix_mus = NULL, mix_sigmas = NULL,
                            mix_skews = NULL, mix_skurts = NULL,
                            mix_fifths = NULL, mix_sixths = NULL) {
  if (is.null(mix_fifths)) {
    e1 <- sum(mix_pis * mix_mus)
    e2 <- sum(mix_pis * (mix_sigmas^2 + mix_mus^2))
    e3 <- sum(mix_pis * (mix_sigmas^3 * mix_skews + 3 * mix_sigmas^2 *
      mix_mus + mix_mus^3))
    e4 <- sum(mix_pis * (mix_sigmas^4 * (mix_skurts + 3) + 4 *
      mix_sigmas^3 * mix_mus * mix_skews + 6 * mix_sigmas^2 * mix_mus^2 +
      mix_mus^4))
    mu3 <- e3 - 3 * e1 * e2 + 2 * e1^3
    mu4 <- e4 - 4 * e1 * e3 + 6 * e1^2 * e2 - 3 * e1^4
    Var <- sum(mix_pis * (mix_sigmas^2 + mix_mus^2)) -
      (sum(mix_pis * mix_mus))^2
    g1 <- mu3/(Var^(3/2))
    g2 <- mu4/(Var^2) - 3
    stcums <- c(e1, sqrt(Var), g1, g2)
    names(stcums) <- c("mean", "sd", "skew", "kurtosis")
    return(stcums)
  }
  e1 <- sum(mix_pis * mix_mus)
  e2 <- sum(mix_pis * (mix_sigmas^2 + mix_mus^2))
  e3 <- sum(mix_pis * (mix_sigmas^3 * mix_skews + 3 * mix_sigmas^2 * mix_mus +
    mix_mus^3))
  e4 <- sum(mix_pis * (mix_sigmas^4 * (mix_skurts + 3) + 4 * mix_sigmas^3 *
    mix_mus * mix_skews + 6 * mix_sigmas^2 * mix_mus^2 + mix_mus^4))
  e5 <- sum(mix_pis * (mix_sigmas^5 * (mix_fifths + 10 * mix_skews) + 5 *
    mix_sigmas ^ 4 * mix_mus * (mix_skurts + 3) + 10 * mix_sigmas ^ 3 *
    mix_mus^2 * mix_skews + 10 * mix_sigmas^2 * mix_mus^3 + mix_mus^5))
  e6 <- sum(mix_pis * (mix_sigmas^6 * (mix_sixths + 15 * mix_skurts + 10 *
    mix_skews^2 + 15) + 6 * mix_sigmas^5 * mix_mus * (mix_fifths + 10 *
    mix_skews) + 15 * mix_sigmas^4 * mix_mus^2 * (mix_skurts + 3) + 20 *
    mix_sigmas^3 * mix_mus^3 * mix_skews + 15 * mix_sigmas^2 * mix_mus^4 +
    mix_mus^6))
  mu3 <- e3 - 3 * e1 * e2 + 2 * e1^3
  mu4 <- e4 - 4 * e1 * e3 + 6 * e1^2 * e2 - 3 * e1^4
  mu5 <- e5 - 5 * e1 * e4 + 10 * e1^2 * e3 - 10 * e1^3 * e2 + 4 * e1^5
  mu6 <- e6 - 6 * e1 * e5 + 15 * e1^2 * e4 - 20 * e1^3 * e3 +
    15 * e1^4 * e2 - 5 * e1^6
  Var <- sum(mix_pis * (mix_sigmas^2 + mix_mus^2)) -
    (sum(mix_pis * mix_mus))^2
  g1 <- mu3/(Var^(3/2))
  g2 <- mu4/(Var^2) - 3
  g3 <- mu5/(Var^(5/2)) - 10 * g1
  g4 <- mu6/(Var^3) - 15 * g2 - 10 * g1^2 - 15
  stcums <- c(e1, sqrt(Var), g1, g2, g3, g4)
  names(stcums) <- c("mean", "sd", "skew", "kurtosis", "fifth", "sixth")
  return(stcums)
}
