#' @title Generation of One Continuous Variable with a Mixture Distribution Using the Power Method Transformation
#'
#' @description This function simulates one continuous mixture variable.  Mixture distributions describe random variables that
#'     are drawn from more than one component distribution.  For a random variable \eqn{Y_{mix}} from a finite continuous mixture
#'     distribution with \eqn{k} components, the probability density function (PDF) can be described by:
#'
#'     \deqn{h_Y(y) = \sum_{i=1}^{k} \pi_i f_{Yi}(y), \sum_{i=1}^{k} \pi_i = 1.}
#'
#'     The \eqn{\pi_i} are mixing parameters which determine the weight of each component distribution \eqn{f_{Yi}(y)} in the overall
#'     probability distribution.  As long as each component has a valid PDF, the overall distribution \eqn{h_Y(y)} has a valid PDF.
#'     The main assumption is statistical independence between the process of randomly selecting the component distribution and the
#'     distributions themselves.  Each component \eqn{Y_i} is generated using either Fleishman's third-order (\code{method} = "Fleishman",
#'     \doi{10.1007/BF02293811}) or Headrick's fifth-order (\code{method} = "Polynomial",
#'     \doi{10.1016/S0167-9473(02)00072-5}) power method transformation (PMT).  It works by matching standardized
#'     cumulants -- the first four (mean, variance, skew, and standardized kurtosis) for Fleishman's method, or the first six (mean,
#'     variance, skew, standardized kurtosis, and standardized fifth and sixth cumulants) for Headrick's method.  The transformation is
#'     expressed as follows:
#'
#'     \deqn{Y = c_0 + c_1 * Z + c_2 * Z^2 + c_3 * Z^3 + c_4 * Z^4 + c_5 * Z^5,  Z \sim N(0,1),}
#'
#'     where \eqn{c_4} and \eqn{c_5} both equal \eqn{0} for Fleishman's method.  The real constants are calculated by \cr
#'     \code{\link[SimMultiCorrData]{find_constants}}.  These components are then transformed to the desired mixture variable using a
#'     random multinomial variable generated based on the mixing probabilities.  There are no parameter input checks in order to decrease
#'     simulation time.  All inputs should be checked prior to simulation with \code{\link[SimCorrMix]{validpar}}.  Summaries for the
#'     simulation results can be obtained with \code{\link[SimCorrMix]{summary_var}}.
#'
#'     Mixture distributions provide a useful way for describing heterogeneity in a population, especially when an outcome is a
#'     composite response from multiple sources.  The vignette \bold{Variable Types} provides more information about simulation of mixture
#'     variables and the required parameters.  The vignette \bold{Expected Cumulants and Correlations for Continuous Mixture Variables}
#'     gives the equations for the expected cumulants of a mixture variable.  In addition, Headrick & Kowalchuk (2007,
#'     \doi{10.1080/10629360600605065}) outlined a general method for comparing a simulated distribution \eqn{Y} to a given theoretical
#'     distribution \eqn{Y^*}.  These steps can be found in the \bold{Continuous Mixture Distributions} vignette.
#'
#' @section Choice of Fleishman's third-order or Headrick's fifth-order method:
#'     Using the fifth-order approximation allows additional control over the fifth and sixth moments of the generated distribution, improving
#'     accuracy.  In addition, the range of feasible standardized kurtosis values, given skew and standardized fifth (\eqn{\gamma_3}) and sixth
#'     (\eqn{\gamma_4}) cumulants, is larger than with Fleishman's method (see \code{\link[SimMultiCorrData]{calc_lower_skurt}}).
#'     For example, the Fleishman method can not be used to generate a non-normal distribution with a ratio of
#'     \eqn{\gamma_3^2/\gamma_4 > 9/14} (see Headrick & Kowalchuk, 2007).  This eliminates the Chi-squared family of distributions, which has
#'     a constant ratio of \eqn{\gamma_3^2/\gamma_4 = 2/3}.  The fifth-order method also generates more distributions with valid PDF's.
#'     However, if the fifth and sixth cumulants are unknown or do not exist, the Fleishman approximation should be used.
#'
#' @section Overview of Simulation Process:
#'     1) A check is performed to see if any distributions are repeated within the parameter inputs, i.e. if the mixture variable
#'     contains 2 components with the same standardized cumulants.  These are noted so that the constants are only calculated once.
#'
#'     2) The constants are calculated for each component variable using \code{\link[SimMultiCorrData]{find_constants}}.  If no
#'     solutions are found that generate a valid power method PDF, the function will return constants that produce an invalid PDF
#'     (or a stop error if no solutions can be found).  Possible solutions include: 1) changing the seed, or 2) using a \code{mix_Six}
#'     list with vectors of sixth cumulant correction values (if \code{method} = "Polynomial").  Errors regarding constant
#'     calculation are the most probable cause of function failure.
#'
#'     3) A matrix \code{X_cont} of dim \code{n x length(mix_pis)} of standard normal variables is generated and singular-value decomposition is done to
#'     remove any correlation.  The \code{constants} are applied to \code{X_cont} to create the component variables \code{Y} with the desired distributions.
#'
#'     4) A random multinomial variable \code{M = rmultinom(n, size = 1, prob = mix_pis)} is generated using \code{\link[stats;Multinom]{rmultinom}}.
#'     The continuous mixture variable \code{Y_mix} is created from the component variables \code{Y} based on this multinomial variable.
#'     That is, if \code{M[i, k_i] = 1}, then \code{Y_mix[i] = Y[i, k_i]}.  A location-scale transformation is done on \code{Y_mix} to give it mean \code{means} and variance \code{vars}.
#'
#' @section Reasons for Function Errors:
#'     1) The most likely cause for function errors is that no solutions to \code{\link[SimMultiCorrData]{fleish}} or
#'     \code{\link[SimMultiCorrData]{poly}} converged when using \code{\link[SimMultiCorrData]{find_constants}}.  If this happens,
#'     the simulation will stop.  It may help to first use \code{\link[SimMultiCorrData]{find_constants}} for each component variable to
#'     determine if a sixth cumulant correction value is needed.  The solutions can be used as starting values (see \code{cstart} below).
#'     If the standardized cumulants are obtained from \code{calc_theory}, the user may need to use rounded values as inputs (i.e.
#'     \code{skews = round(skews, 8)}).  For example, in order to ensure that skew is exactly 0 for symmetric distributions.
#'
#'     2) The kurtosis may be outside the region of possible values.  There is an associated lower boundary for kurtosis associated
#'     with a given skew (for Fleishman's method) or skew and fifth and sixth cumulants (for Headrick's method).  Use
#'     \code{\link[SimMultiCorrData]{calc_lower_skurt}} to determine the boundary for a given set of cumulants.
#'
#' @param n the sample size (i.e. the length of the simulated variable; default = 10000)
#' @param method the method used to generate the component variables.  "Fleishman" uses Fleishman's third-order polynomial transformation
#'     and "Polynomial" uses Headrick's fifth-order transformation.
#' @param means mean for the mixture variable (default = 0)
#' @param vars variance for the mixture variable (default = 1)
#' @param mix_pis a vector of mixing probabilities that sum to 1 for the component distributions
#' @param mix_mus a vector of means for the component distributions
#' @param mix_sigmas a vector of standard deviations for the component distributions
#' @param mix_skews a vector of skew values for the component distributions
#' @param mix_skurts a vector of standardized kurtoses for the component distributions
#' @param mix_fifths a vector of standardized fifth cumulants for the component distributions; keep NULL if using \code{method} = "Fleishman"
#'     to generate continuous variables
#' @param mix_sixths a vector of standardized sixth cumulants for the component distributions; keep NULL if using \code{method} = "Fleishman"
#'     to generate continuous variables
#' @param mix_Six a list of vectors of sixth cumulant correction values for the component distributions of \eqn{Y_{mix}};
#'     use \code{NULL} if no correction is desired for a given component; if no correction is desired for any component keep as
#'     \code{mix_Six = list()} (not necessary for \code{method} = "Fleishman")
#' @param seed the seed value for random number generation (default = 1234)
#' @param cstart a list of length equal to the total number of mixture components containing initial values for root-solving
#'     algorithm used in \code{\link[SimMultiCorrData]{find_constants}}.  If user specified, each list element must be input as a matrix.
#'     For \code{method} = "Fleishman", each should have 3 columns for \eqn{c_1, c_2, c_3};
#'     for \code{method} = "Polynomial", each should have 5 columns for \eqn{c_1, c_2, c_3, c_4, c_5}.  If no starting values are specified for
#'     a given component, that list element should be \code{NULL}.
#' @param quiet if FALSE prints total simulation time
#' @import SimMultiCorrData
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @import BB
#' @import nleqslv
#' @export
#' @keywords simulation continuous mixture Fleishman Headrick
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimCorrMix]{validpar}}, \code{\link[SimCorrMix]{summary_var}}
#' @return A list with the following components:
#' @return \code{constants} a data.frame of the constants
#' @return \code{Y_comp} a data.frame of the components of the mixture variable
#' @return \code{Y_mix} a data.frame of the generated mixture variable
#' @return \code{sixth_correction} the sixth cumulant correction values for \code{Y_comp}
#' @return \code{valid.pdf} "TRUE" if constants generate a valid PDF, else "FALSE"
#' @return \code{Time} the total simulation time in minutes
#' @references
#' Davenport JW, Bezder JC, & Hathaway RJ (1988). Parameter Estimation for Finite Mixture Distributions.
#'     Computers & Mathematics with Applications, 15(10):819-28.
#'
#' Everitt BS (1996). An Introduction to Finite Mixture Distributions. Statistical Methods in Medical Research, 5(2):107-127. \doi{10.1177/096228029600500202}.
#'
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43:521-532. \doi{10.1007/BF02293811}.
#'
#' Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#'     Non-normal Distributions. Computational Statistics & Data Analysis, 40(4):685-711. \doi{10.1016/S0167-9473(02)00072-5}.
#'     (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3(1):65-71. \doi{10.22237/jmasm/1083370080}.
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77:229-249. \doi{10.1080/10629360600605065}.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64:25-35. \doi{10.1007/BF02294317}.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3):1 - 17. \cr \doi{10.18637/jss.v019.i03}.
#'
#' Pearson, RK. 2011. "Exploring Data in Engineering, the Sciences, and Medicine." In. New York: Oxford University Press.
#'
#' @examples \dontrun{
#' # Mixture of Beta(6, 3), Beta(4, 1.5), and Beta(10, 20)
#' Stcum1 <- calc_theory("Beta", c(6, 3))
#' Stcum2 <- calc_theory("Beta", c(4, 1.5))
#' Stcum3 <- calc_theory("Beta", c(10, 20))
#' mix_pis <- c(0.5, 0.2, 0.3)
#' mix_mus <- c(Stcum1[1], Stcum2[1], Stcum3[1])
#' mix_sigmas <- c(Stcum1[2], Stcum2[2], Stcum3[2])
#' mix_skews <- c(Stcum1[3], Stcum2[3], Stcum3[3])
#' mix_skurts <- c(Stcum1[4], Stcum2[4], Stcum3[4])
#' mix_fifths <- c(Stcum1[5], Stcum2[5], Stcum3[5])
#' mix_sixths <- c(Stcum1[6], Stcum2[6], Stcum3[6])
#' mix_Six <- list(seq(0.01, 10, 0.01), c(0.01, 0.02, 0.03),
#'   seq(0.01, 10, 0.01))
#' Bstcum <- calc_mixmoments(mix_pis, mix_mus, mix_sigmas, mix_skews,
#'   mix_skurts, mix_fifths, mix_sixths)
#' Bmix <- contmixvar1(n = 10000, "Polynomial", Bstcum[1], Bstcum[2]^2,
#'   mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, mix_fifths,
#'   mix_sixths, mix_Six)
#' Bsum <- summary_var(Y_comp = Bmix$Y_comp, Y_mix = Bmix$Y_mix, means = means,
#'   vars = vars, mix_pis = mix_pis, mix_mus = mix_mus,
#'   mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts,
#'   mix_fifths = mix_fifths, mix_sixths = mix_sixths)
#' }
contmixvar1 <- function(n = 10000, method = c("Fleishman", "Polynomial"),
                        means = 0, vars = 1, mix_pis = NULL, mix_mus = NULL,
                        mix_sigmas = NULL, mix_skews =  NULL,
                        mix_skurts =  NULL, mix_fifths =  NULL,
                        mix_sixths =  NULL, mix_Six = list(), seed = 1234,
                        cstart = list(), quiet = FALSE) {
  start.time <- Sys.time()
  csame.dist <- NULL
  for (i in 2:length(mix_skews)) {
    if (mix_skews[i] %in% mix_skews[1:(i - 1)]) {
      csame <- which(mix_skews[1:(i - 1)] == mix_skews[i])
      for (j in 1:length(csame)) {
        if (method == "Polynomial") {
          if ((mix_skurts[i] == mix_skurts[csame[j]]) &
              (mix_fifths[i] == mix_fifths[csame[j]]) &
              (mix_sixths[i] == mix_sixths[csame[j]])) {
            csame.dist <- rbind(csame.dist, c(csame[j], i))
            break
          }
        }
        if (method == "Fleishman") {
          if (mix_skurts[i] == mix_skurts[csame[j]]) {
            csame.dist <- rbind(csame.dist, c(csame[j], i))
            break
          }
        }
      }
    }
  }
  SixCorr <- numeric(length(mix_pis))
  Valid.PDF <- numeric(length(mix_pis))
  if (method == "Fleishman") {
    constants <- matrix(NA, nrow = length(mix_pis), ncol = 4)
    colnames(constants) <- c("c0", "c1", "c2", "c3")
  }
  if (method == "Polynomial") {
    constants <- matrix(NA, nrow = length(mix_pis), ncol = 6)
    colnames(constants) <- c("c0", "c1", "c2", "c3", "c4", "c5")
  }
  for (i in 1:length(mix_pis)) {
    if (!is.null(csame.dist)) {
      rind <- which(csame.dist[, 2] == i)
      if (length(rind) > 0) {
        constants[i, ] <- constants[csame.dist[rind, 1], ]
        SixCorr[i] <- SixCorr[csame.dist[rind, 1]]
        Valid.PDF[i] <- Valid.PDF[csame.dist[rind, 1]]
      }
    }
    if (sum(is.na(constants[i, ])) > 0) {
      if (length(mix_Six) == 0) Six2 <- NULL else
        Six2 <- mix_Six[[i]]
      if (length(cstart) == 0) cstart2 <- NULL else
        cstart2 <- cstart[[i]]
      cons <-
        suppressWarnings(find_constants(method = method, skews = mix_skews[i],
          skurts = mix_skurts[i], fifths = mix_fifths[i],
          sixths = mix_sixths[i], Six = Six2, cstart = cstart2, n = 25,
          seed = seed))
      if (length(cons) == 1 | is.null(cons)) {
        stop(paste("Constants can not be found for component ", i,
                   ".", sep = ""))
      }
      con_solution <- cons$constants
      SixCorr[i] <- ifelse(is.null(cons$SixCorr1), NA, cons$SixCorr1)
      Valid.PDF[i] <- cons$valid
      constants[i, ] <- con_solution
    }
  }
  set.seed(seed)
  X_cont <- matrix(rnorm(length(mix_pis) * n), n)
  X_cont <- scale(X_cont, TRUE, FALSE)
  X_cont <- X_cont %*% svd(X_cont, nu = 0)$v
  X_cont <- scale(X_cont, FALSE, TRUE)
  Y <- matrix(1, nrow = n, ncol = length(mix_pis))
  Yb <- matrix(1, nrow = n, ncol = length(mix_pis))
  for (i in 1:length(mix_pis)) {
    if (method == "Fleishman") {
      Y[, i] <- constants[i, 1] + constants[i, 2] * X_cont[, i] +
        constants[i, 3] * X_cont[, i]^2 + constants[i, 4] * X_cont[, i]^3
    }
    if (method == "Polynomial") {
      Y[, i] <- constants[i, 1] + constants[i, 2] * X_cont[, i] +
        constants[i, 3] * X_cont[, i]^2 + constants[i, 4] * X_cont[, i]^3 +
        constants[i, 5] * X_cont[, i]^4 + constants[i, 6] * X_cont[, i]^5
    }
    Yb[, i] <- mix_mus[i] + mix_sigmas[i] * Y[, i]
  }
  set.seed(seed)
  M <- rmultinom(n, size = 1, prob = mix_pis)
  Y_mix <- apply(t(M) * Yb, 1, sum)
  Y_mix <- scale(Y_mix)
  Y_mix <- matrix(means + sqrt(vars) * Y_mix, n, 1)
  stop.time <- Sys.time()
  Time <- round(difftime(stop.time, start.time, units = "min"), 3)
  if (quiet == FALSE) cat("Total Simulation time:", Time, "minutes \n")
  result <- list(constants = as.data.frame(constants),
    Y_comp = Yb, Y_mix = Y_mix, sixth_correction = SixCorr,
    valid.pdf = Valid.PDF, Time = Time)
  result
}
