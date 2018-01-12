#' @title Parameter Check for Simulation or Correlation Validation Functions
#'
#' @description This function checks the parameter inputs to the simulation functions \code{\link[SimCorrMix]{contmixvar1}},
#'     \code{\link[SimCorrMix]{corrvar}}, and \code{\link[SimCorrMix]{corrvar2}} and to the correlation validation functions
#'     \code{\link[SimCorrMix]{validcorr}} and \code{\link[SimCorrMix]{validcorr2}}.  It should be used prior to execution of these
#'     functions to ensure all inputs are of the correct format.  Those functions do not contain parameter checks in order to decrease
#'     simulation time.  This would be important if the user is running several simulation repetitions so that the inputs only have to
#'     be checked once.  Note that the inputs do not include all of the inputs to the simulation functions.  See the appropriate function
#'     documentation for more details about parameter inputs.
#'
#' @param k_cat the number of ordinal (r >= 2 categories) variables (default = 0)
#' @param k_cont the number of continuous non-mixture variables (default = 0)
#' @param k_mix the number of continuous mixture variables (default = 0)
#' @param k_pois the number of regular Poisson and zero-inflated Poisson variables (default = 0)
#' @param k_nb the number of regular Negative Binomial and zero-inflated Negative Binomial variables (default = 0)
#' @param method the method used to generate the \code{k_cont} non-mixture and \code{k_mix} mixture continuous variables.
#'     "Fleishman" uses Fleishman's third-order polynomial transformation and "Polynomial" uses Headrick's fifth-order transformation.
#' @param means a vector of means for the \code{k_cont} non-mixture and \code{k_mix} mixture continuous variables
#'     (i.e. \code{rep(0, (k_cont + k_mix))})
#' @param vars a vector of variances for the \code{k_cont} non-mixture and \code{k_mix} mixture continuous variables
#'     (i.e. \code{rep(1, (k_cont + k_mix))})
#' @param skews a vector of skewness values for the \code{k_cont} non-mixture continuous variables
#' @param skurts a vector of standardized kurtoses (kurtosis - 3, so that normal variables have a value of 0)
#'     for the \code{k_cont} non-mixture continuous variables
#' @param fifths a vector of standardized fifth cumulants for the \code{k_cont} non-mixture continuous variables
#'     (not necessary for \code{method} = "Fleishman")
#' @param sixths a vector of standardized sixth cumulants for the \code{k_cont} non-mixture continuous variables
#'     (not necessary for \code{method} = "Fleishman")
#' @param Six a list of vectors of sixth cumulant correction values for the \code{k_cont} non-mixture continuous variables
#'     if no valid PDF constants are found, \cr ex: \code{Six = list(seq(0.01, 2, 0.01), seq(1, 10, 0.5))};
#'     if no correction is desired for variable \eqn{Y_{cont_i}}, set set the i-th list component equal to \code{NULL};
#'     if no correction is desired for any of the \eqn{Y_{cont}} keep as \code{Six = list()}
#'     (not necessary for \code{method} = "Fleishman")
#' @param mix_pis a vector if using \code{\link[SimCorrMix]{contmixvar1}} or a list of length \code{k_mix} with i-th component a vector of mixing probabilities that sum to 1 for component distributions of \eqn{Y_{mix_i}}
#' @param mix_mus a vector if using \code{\link[SimCorrMix]{contmixvar1}} or a list of length \code{k_mix} with i-th component a vector of means for component distributions of \eqn{Y_{mix_i}}
#' @param mix_sigmas a vector if using \code{\link[SimCorrMix]{contmixvar1}} or a list of length \code{k_mix} with i-th component a vector of standard deviations for component distributions of \eqn{Y_{mix_i}}
#' @param mix_skews a vector if using \code{\link[SimCorrMix]{contmixvar1}} or a list of length \code{k_mix} with i-th component a vector of skew values for component distributions of \eqn{Y_{mix_i}}
#' @param mix_skurts a vector if using \code{\link[SimCorrMix]{contmixvar1}} or a list of length \code{k_mix} with i-th component a vector of standardized kurtoses for component distributions of \eqn{Y_{mix_i}}
#' @param mix_fifths a vector if using \code{\link[SimCorrMix]{contmixvar1}} or a list of length \code{k_mix} with i-th component a vector of standardized fifth cumulants for component distributions of \eqn{Y_{mix_i}}
#'     (not necessary for \code{method} = "Fleishman")
#' @param mix_sixths a vector if using \code{\link[SimCorrMix]{contmixvar1}} or a list of length \code{k_mix} with i-th component a vector of standardized sixth cumulants for component distributions of \eqn{Y_{mix_i}}
#'     (not necessary for \code{method} = "Fleishman")
#' @param mix_Six if using \code{\link[SimCorrMix]{contmixvar1}}, a list of vectors of sixth cumulant corrections for the components of the
#'     continuous mixture variable; else a list of length \code{k_mix} with i-th component a list of vectors of sixth cumulant correction values
#'     for component distributions of \eqn{Y_{mix_i}}; use \code{NULL} if no correction is desired for a given component or
#'     mixture variable; if no correction is desired for any of the \eqn{Y_{mix}} keep as \code{mix_Six = list()}
#'     (not necessary for \code{method} = "Fleishman")
#' @param marginal a list of length equal to \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1; default = list());
#'     for binary variables, these should be input the same as for ordinal variables with more than 2 categories (i.e. the user-specified
#'     probability is the probability of the 1st category, which has the smaller support value)
#' @param support a list of length equal to \code{k_cat}; the i-th element is a vector containing the r ordered support values;
#'     if not provided (i.e. \code{support = list()}), the default is for the i-th element to be the vector 1, ..., r
#' @param lam a vector of lambda (mean > 0) constants for the Poisson variables (see \code{\link[stats;Poisson]{dpois}}); the order should be
#'     1st regular Poisson variables, 2nd zero-inflated Poisson variables
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
#' @param rho the target correlation matrix which must be ordered
#'     \emph{1st ordinal, 2nd continuous non-mixture, 3rd components of continuous mixtures, 4th regular Poisson, 5th zero-inflated Poisson,
#'     6th regular NB, 7th zero-inflated NB}; note that \code{rho} is specified in terms of the components of \eqn{Y_{mix}}
#' @param Sigma an intermediate correlation matrix to use if the user wants to provide one, else it is calculated within by
#'     \code{\link[SimCorrMix]{intercorr}}
#' @param cstart a list of length equal to \code{k_cont} + the total number of mixture components containing initial values for root-solving
#'     algorithm used in \code{\link[SimMultiCorrData]{find_constants}}.  If user specified, each list element must be input as a matrix.
#'     For \code{method} = "Fleishman", each should have 3 columns for \eqn{c1, c2, c3};
#'     for \code{method} = "Polynomial", each should have 5 columns for \eqn{c1, c2, c3, c4, c5}.  If no starting values are specified for
#'     a given component, that list element should be \code{NULL}.
#' @param quiet if FALSE prints messages, if TRUE suppresses message printing
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @export
#' @return TRUE if all inputs are correct, else it will stop with a correction message
#' @keywords ParameterCheck
#' @seealso \code{\link[SimCorrMix]{contmixvar1}}, \code{\link[SimCorrMix]{corrvar}}, \code{\link[SimCorrMix]{corrvar2}},
#'     \code{\link[SimCorrMix]{validcorr}}, \code{\link[SimCorrMix]{validcorr2}}
#' @examples \dontrun{
#' # 2 continuous mixture, 1 binary, 1 zero-inflated Poisson, and
#' # 1 zero-inflated NB variable
#'
#' # Mixture variables: Normal mixture with 2 components;
#' # mixture of Logistic(0, 1), Chisq(4), Beta(4, 1.5)
#' # Find cumulants of components of 2nd mixture variable
#' L <- calc_theory("Logistic", c(0, 1))
#' C <- calc_theory("Chisq", 4)
#' B <- calc_theory("Beta", c(4, 1.5))
#'
#' skews <- skurts <- fifths <- sixths <- NULL
#' Six <- list()
#' mix_pis <- list(c(0.4, 0.6), c(0.3, 0.2, 0.5))
#' mix_mus <- list(c(-2, 2), c(L[1], C[1], B[1]))
#' mix_sigmas <- list(c(1, 1), c(L[2], C[2], B[2]))
#' mix_skews <- list(rep(0, 2), c(L[3], C[3], B[3]))
#' mix_skurts <- list(rep(0, 2), c(L[4], C[4], B[4]))
#' mix_fifths <- list(rep(0, 2), c(L[5], C[5], B[5]))
#' mix_sixths <- list(rep(0, 2), c(L[6], C[6], B[6]))
#' mix_Six <- list(list(NULL, NULL), list(1.75, NULL, 0.03))
#' Nstcum <- calc_mixmoments(mix_pis[[1]], mix_mus[[1]], mix_sigmas[[1]],
#'   mix_skews[[1]], mix_skurts[[1]], mix_fifths[[1]], mix_sixths[[1]])
#' Mstcum <- calc_mixmoments(mix_pis[[2]], mix_mus[[2]], mix_sigmas[[2]],
#'   mix_skews[[2]], mix_skurts[[2]], mix_fifths[[2]], mix_sixths[[2]])
#' means <- c(Nstcum[1], Mstcum[1])
#' vars <- c(Nstcum[2]^2, Mstcum[2]^2)
#'
#' marginal <- list(0.3)
#' support <- list(c(0, 1))
#' lam <- 0.5
#' p_zip <- 0.1
#' size <- 2
#' prob <- 0.75
#' p_zinb <- 0.2
#'
#' k_cat <- k_pois <- k_nb <- 1
#' k_cont <- 0
#' k_mix <- 2
#' Rey <- matrix(0.39, 8, 8)
#' diag(Rey) <- 1
#' rownames(Rey) <- colnames(Rey) <- c("O1", "M1_1", "M1_2", "M2_1", "M2_2",
#'   "M2_3", "P1", "NB1")
#'
#' # set correlation between components of the same mixture variable to 0
#' Rey["M1_1", "M1_2"] <- Rey["M1_2", "M1_1"] <- 0
#' Rey["M2_1", "M2_2"] <- Rey["M2_2", "M2_1"] <- Rey["M2_1", "M2_3"] <- 0
#' Rey["M2_3", "M2_1"] <- Rey["M2_2", "M2_3"] <- Rey["M2_3", "M2_2"] <- 0
#'
#' # use before contmixvar1 with 1st mixture variable:
#' # change mix_pis to not sum to 1
#'
#' check1 <- validpar(k_mix = 1, method = "Polynomial", means = Nstcum[1],
#'   vars = Nstcum[2]^2, mix_pis = C(0.4, 0.5), mix_mus = mix_mus[[1]],
#'   mix_sigmas = mix_sigmas[[1]], mix_skews = mix_skews[[1]],
#'   mix_skurts = mix_skurts[[1]], mix_fifths = mix_fifths[[1]],
#'   mix_sixths = mix_sixths[[1]])
#'
#' # use before validcorr: should return TRUE
#'
#' check2 <- validpar(k_cat, k_cont, k_mix, k_pois, k_nb, "Polynomial", means,
#'   vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
#'   mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support,
#'   lam, p_zip, size, prob, mu = NULL, p_zinb, rho = Rey)
#'
#' }
#'
validpar <- function(k_cat = 0, k_cont = 0, k_mix = 0, k_pois = 0,
                     k_nb = 0, method = c("Fleishman", "Polynomial"),
                     means =  NULL, vars =  NULL, skews =  NULL,
                     skurts =  NULL, fifths =  NULL, sixths =  NULL,
                     Six = list(), mix_pis = list(), mix_mus = list(),
                     mix_sigmas = list(), mix_skews =  list(),
                     mix_skurts =  list(), mix_fifths =  list(),
                     mix_sixths =  list(), mix_Six = list(),
                     marginal = list(), support = list(), lam  =  NULL,
                     p_zip = 0, size = NULL, prob = NULL, mu = NULL,
                     p_zinb = 0, pois_eps = 0.0001, nb_eps = 0.0001,
                     rho = NULL, Sigma = NULL, cstart = list(),
                     quiet = FALSE) {
  if (k_cat > 0) {
    if (k_cat != length(marginal))
      stop("Length of marginal does not match the number of ordinal
           variables.")
    if (!all(unlist(lapply(marginal,
        function(x) (sort(x) == x & min(x) > 0 & max(x) < 1)))))
      stop("Error in given marginal distributions.")
    if (length(support) > 0 & (k_cat != length(support)))
      stop("Length of support does not match the number of ordinal
           variables.")
  }
  k_comp <- 0
  if ((k_cont + k_mix) > 0) {
    if (length(method) != 1)
      stop("Choose a PMT method for the continuous variables.")
    if (length(means) != (k_cont + k_mix) | length(vars) != (k_cont + k_mix))
      stop("Length of means and vars should be k_cont + k_mix.")
    if (k_cont > 0) {
      if (length(skews) != k_cont | length(skurts) != k_cont)
        stop("Parameters for continuous non-mixture distributions should be
             vectors of length k_cont.")
      if (method == "Polynomial") {
        if (length(fifths) != k_cont | length(sixths) != k_cont |
            (length(Six) != 0 & length(Six) != k_cont))
          stop("Parameters for continuous non-mixture distributions should be
               vectors of length k_cont. Six should be either list() or
               a list of length k_cont.")
      }
    }
    if (k_mix > 0) {
      if (class(mix_pis) == "list") {
        k_comp <- lengths(mix_pis)
        if (length(mix_pis) != k_mix | length(mix_mus) != k_mix |
            length(mix_sigmas) != k_mix | length(mix_skews) != k_mix |
            length(mix_skurts) != k_mix)
          stop("Parameters for continuous mixture distributions should be
               lists of length equal to k_mix.")
        if (!all(lengths(mix_mus) %in% k_comp) |
            !all(lengths(mix_sigmas) %in% k_comp) |
            !all(lengths(mix_skews) %in% k_comp) |
            !all(lengths(mix_skurts) %in% k_comp))
          stop("Components of mixture parameter lists should be the same
               length as components of mix_pis.")
        if (method == "Polynomial") {
          if (length(mix_fifths) != k_mix | length(mix_sixths) != k_mix |
              (length(mix_Six) != 0 & length(mix_Six) != k_mix))
            stop("Parameters for continuous mixture distributions should be
                 lists of length equal to k_mix. mix_Six should be either
                 list() or a list of length k_mix.")
          if (!all(lengths(mix_fifths) %in% k_comp) |
              !all(lengths(mix_sixths) %in% k_comp))
            stop("Components of mixture parameter lists should be the same
                 length as components of mix_pis.")
        }
        if (all.equal(unlist(lapply(mix_pis, sum)), rep(1, k_mix)) == FALSE)
          stop("Mixing parameters should sum to 1 for each variable.")
      }
      if (class(mix_pis) == "numeric") {
        if (k_mix > 1 | k_cat > 0 | k_pois > 0 | k_nb > 0)
          stop("Mixture parameters should be lists of length equal to k_mix.")
        if (sum(mix_pis) != 1)
          stop("Mixing parameters should sum to 1.")
        k_comp <- length(mix_pis)
        if (length(mix_mus) != k_comp | length(mix_sigmas) != k_comp |
            length(mix_skews) != k_comp | length(mix_skurts) != k_comp)
          stop("Mixture parameters should have same length as mix_pis.")
        if (method == "Polynomial") {
          if (length(mix_fifths) != k_comp | length(mix_sixths) != k_comp |
              !(length(mix_Six) %in% c(0, k_comp)))
            stop("Mixture parameters should have same length as mix_pis.")
        }
      }
    }
    if (length(cstart) != 0) {
      if (length(cstart) != (k_cont + k_comp))
        stop("Length of cstart should be k_cont + number of mixture
             components.")
      k.c <- unlist(lapply(cstart, ncol))
      if (method == "Fleishman") {
        if (!all.equal(k.c, rep(3, length(cstart))))
          stop("Dimension of cstart matrices should be number of starting
               values by 3")
      } else {
        if (!all.equal(k.c, rep(5, length(cstart))))
          stop("Dimension of cstart matrices should be number of starting
               values by 5")
      }
    }
  }
  if (k_pois > 0) {
    if (k_pois != length(lam))
      stop("Length of lam does not match the number of Poisson variables.")
    if (sum(lam < 0) > 0)
      stop("Lambda values cannot be negative.")
    if (length(p_zip) < k_pois & quiet == FALSE)
      message("Default of p_zip = 0 will be used for Poisson variables.")
    if (length(pois_eps) < k_pois & quiet == FALSE)
      message("Default of pois_eps = 0.0001 will be used for Poisson variables
              if using correlation method 2.")
  }
  if (k_nb > 0) {
    if (k_nb != length(size))
      stop("Length of size does not match the number of NB variables.")
    if (length(prob) > 0 & length(mu) > 0)
      stop("Either give success probabilities or means for NB variables.")
    if (length(prob) > 0 & k_nb != length(prob))
      stop("Length of prob does not match the number of NB variables.")
    if (length(mu) > 0 & k_nb != length(mu))
      stop("Length of mu does not match the number of NB variables.")
    if (length(p_zinb) < k_nb & quiet == FALSE)
      message("Default of p_zinb = 0 will be used for NB variables.")
    if (length(nb_eps) < k_nb & quiet == FALSE)
      message("Default of nb_eps = 0.0001 will be used for NB variables
              if using correlation method 2.")
  }
  k <- k_cat + k_cont + sum(k_comp) + k_pois + k_nb
  if (!is.null(Sigma)) {
    if (ncol(Sigma) != k | nrow(Sigma) != k)
      stop("Sigma matrix is not of the right dimension.")
    if (!isSymmetric(Sigma) | !all(diag(Sigma) == 1))
      stop("Sigma matrix not valid! Check symmetry and diagonal values.")
    if (min(eigen(Sigma, symmetric = TRUE)$values) < 0 & quiet == FALSE)
      message("Sigma matrix is not positive-definite.  Nearest
              positive-definite matrix will be used if use.nearPD = TRUE.")
  }
  if (!is.null(rho)) {
    if (ncol(rho) != k | nrow(rho) != k)
      stop("Correlation matrix is not of the right dimension.")
    if (!isSymmetric(rho) | !all(diag(rho) == 1))
      stop("Correlation matrix not valid! Check symmetry and diagonal
           values.")
    if (min(eigen(rho, symmetric = TRUE)$values) < 0 & quiet == FALSE)
      message("Target correlation matrix is not positive definite.")
  }
  return(TRUE)
}
