#' @title Summary of Simulated Variables
#'
#' @description This function summarizes the results of \code{\link[SimCorrMix]{contmixvar1}}, \code{\link[SimCorrMix]{corrvar}}, or
#'     \code{\link[SimCorrMix]{corrvar2}}.  The inputs are either the simulated variables or inputs for those functions.  See their
#'     documentation for more information.  If summarizing result from \code{\link[SimCorrMix]{contmixvar1}}, mixture parameters may be
#'     entered as vectors instead of lists.
#'
#' @param Y_cat a matrix of ordinal variables
#' @param Y_cont a matrix of continuous non-mixture variables
#' @param Y_comp a matrix of components of continuous mixture variables
#' @param Y_mix a matrix of continuous mixture variables
#' @param Y_pois a matrix of Poisson variables
#' @param Y_nb a matrix of Negative Binomial variables
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
#' @param mix_pis a list of length \code{k_mix} with i-th component a vector of mixing probabilities that sum to 1 for component distributions of \eqn{Y_{mix_i}}
#' @param mix_mus a list of length \code{k_mix} with i-th component a vector of means for component distributions of \eqn{Y_{mix_i}}
#' @param mix_sigmas a list of length \code{k_mix} with i-th component a vector of standard deviations for component distributions of \eqn{Y_{mix_i}}
#' @param mix_skews a list of length \code{k_mix} with i-th component a vector of skew values for component distributions of \eqn{Y_{mix_i}}
#' @param mix_skurts a list of length \code{k_mix} with i-th component a vector of standardized kurtoses for component distributions of \eqn{Y_{mix_i}}
#' @param mix_fifths a list of length \code{k_mix} with i-th component a vector of standardized fifth cumulants for component distributions of \eqn{Y_{mix_i}}
#'     (not necessary for \code{method} = "Fleishman")
#' @param mix_sixths a list of length \code{k_mix} with i-th component a vector of standardized sixth cumulants for component distributions of \eqn{Y_{mix_i}}
#'     (not necessary for \code{method} = "Fleishman")
#' @param marginal a list of length equal to \code{k_cat}; the i-th element is a vector of the cumulative
#'     probabilities defining the marginal distribution of the i-th variable;
#'     if the variable can take r values, the vector will contain r - 1 probabilities (the r-th is assumed to be 1);
#'     for binary variables, these should be input the same as for ordinal variables with more than 2 categories (i.e. the user-specified
#'     probability is the probability of the 1st category, which has the smaller support value)
#' @param lam a vector of lambda (mean > 0) constants for the Poisson variables (see \code{\link[stats]{dpois}}); the order should be
#'     1st regular Poisson variables, 2nd zero-inflated Poisson variables
#' @param p_zip a vector of probabilities of structural zeros (not including zeros from the Poisson distribution) for the
#'     zero-inflated Poisson variables (see \code{\link[VGAM]{dzipois}}); if \code{p_zip} = 0, \eqn{Y_{pois}} has a regular Poisson
#'     distribution; if \code{p_zip} is in (0, 1), \eqn{Y_{pois}} has a zero-inflated Poisson distribution;
#'     if \code{p_zip} is in \code{(-(exp(lam) - 1)^(-1), 0)}, \eqn{Y_{pois}} has a zero-deflated Poisson distribution and \code{p_zip}
#'     is not a probability; if \code{p_zip = -(exp(lam) - 1)^(-1)}, \eqn{Y_{pois}} has a positive-Poisson distribution
#'     (see \code{\link[VGAM]{dpospois}}); if \code{length(p_zip) < length(lam)}, the missing values are set to 0 (and ordered 1st)
#' @param size a vector of size parameters for the Negative Binomial variables (see \code{\link[stats]{dnbinom}}); the order should be
#'     1st regular NB variables, 2nd zero-inflated NB variables
#' @param prob a vector of success probability parameters for the NB variables; order the same as in \code{size}
#' @param mu a vector of mean parameters for the NB variables (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables,
#'     not a mixture; default = NULL); order the same as in \code{size}; for zero-inflated NB this refers to
#'     the mean of the NB distribution (see \code{\link[VGAM]{dzinegbin}})
#' @param p_zinb a vector of probabilities of structural zeros (not including zeros from the NB distribution) for the zero-inflated NB variables
#'     (see \code{\link[VGAM]{dzinegbin}}); if \code{p_zinb} = 0, \eqn{Y_{nb}} has a regular NB distribution;
#'     if \code{p_zinb} is in \code{(-prob^size/(1 - prob^size),} \code{0)}, \eqn{Y_{nb}} has a zero-deflated NB distribution and \code{p_zinb}
#'     is not a probability; if \code{p_zinb = -prob^size/(1 - prob^size)}, \eqn{Y_{nb}} has a positive-NB distribution (see
#'     \code{\link[VGAM]{dposnegbin}}); if \code{length(p_zinb) < length(size)}, the missing values are set to 0 (and ordered 1st)
#' @param rho the target correlation matrix which must be ordered
#'     \emph{1st ordinal, 2nd continuous non-mixture, 3rd components of continuous mixtures, 4th regular Poisson, 5th zero-inflated Poisson,
#'     6th regular NB, 7th zero-inflated NB}; note that \code{rho} is specified in terms of the components of \eqn{Y_{mix}}
#' @importFrom psych describe
#' @import SimMultiCorrData
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @export
#' @keywords summary
#' @seealso \code{\link[SimCorrMix]{contmixvar1}}, \code{\link[SimCorrMix]{corrvar}}, \code{\link[SimCorrMix]{corrvar2}}
#' @return A list whose components vary based on the type of simulated variables.
#' @return If \bold{ordinal variables} are produced:
#'
#'     \code{ord_sum} a list, where the i-th element contains a data.frame with target and simulated cumulative probabilities for ordinal variable Y_i
#' @return If \bold{continuous variables} are produced:
#'
#'     \code{cont_sum} a data.frame summarizing \code{Y_cont} and \code{Y_comp},
#'
#'     \code{target_sum} a data.frame with the target distributions for \code{Y_cont} and \code{Y_comp},
#'
#'     \code{mix_sum} a data.frame summarizing \code{Y_mix},
#'
#'     \code{target_mix} a data.frame with the target distributions for \code{Y_mix},
#'
#' @return If \bold{Poisson variables} are produced:
#'
#'     \code{pois_sum} a data.frame summarizing \code{Y_pois}
#' @return If \bold{Negative Binomial variables} are produced:
#'
#'     \code{nb_sum} a data.frame summarizing \code{Y_nb}
#' @return Additionally, the following elements:
#'
#'     \code{rho_calc} the final correlation matrix for \code{Y_cat}, \code{Y_cont}, \code{Y_comp}, \code{Y_pois}, and \code{Y_nb}
#'
#'     \code{rho_mix} the final correlation matrix for \code{Y_cat}, \code{Y_cont}, \code{Y_mix}, \code{Y_pois}, and \code{Y_nb}
#'
#'     \code{maxerr} the maximum final correlation error of \code{rho_calc} from the target \code{rho}.
#'
#' @references See references for \code{\link[SimCorrMix]{SimCorrMix}}.
#'
#' @examples \dontrun{
#'
#' # 2 continuous mixture, 1 binary, 1 zero-inflated Poisson, and
#' # 1 zero-inflated NB variable
#' n <- 10000
#' seed <- 1234
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
#' # check parameter inputs
#' validpar(k_cat, k_cont, k_mix, k_pois, k_nb, "Polynomial", means,
#'   vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
#'   mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support,
#'   lam, p_zip, size, prob, mu = NULL, p_zinb, rho = Rey)
#'
#' # check to make sure Rey is within the feasible correlation boundaries
#' validcorr(n, k_cat, k_cont, k_mix, k_pois, k_nb, "Polynomial", means,
#'   vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
#'   mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
#'   lam, p_zip, size, prob, mu = NULL, p_zinb, Rey, seed)
#'
#' # simulate without the error loop
#' Sim1 <- corrvar(n, k_cat, k_cont, k_mix, k_pois, k_nb, "Polynomial", means,
#'   vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
#'   mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support,
#'   lam, p_zip, size, prob, mu = NULL, p_zinb, Rey, seed, epsilon = 0.01)
#'
#' Summ1 <- summary_var(Sim1$Y_cat, Y_cont = NULL, Sim1$Y_comp, Sim1$Y_mix,
#'   Sim1$Y_pois, Sim1$Y_nb, means, vars, skews, skurts, fifths, sixths,
#'   mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, mix_fifths,
#'   mix_sixths, marginal, lam, p_zip, size, prob, mu = NULL, p_zinb, Rey)
#'
#' Sim1_error <- abs(Rey - Summ1$rho_calc)
#' summary(as.numeric(Sim1_error))
#'
#' }
#'
#'
summary_var <- function(Y_cat = NULL, Y_cont = NULL, Y_comp = NULL,
                        Y_mix = NULL, Y_pois = NULL, Y_nb = NULL,
                        means =  NULL, vars =  NULL, skews =  NULL,
                        skurts =  NULL, fifths =  NULL, sixths =  NULL,
                        mix_pis = list(), mix_mus = list(),
                        mix_sigmas = list(), mix_skews =  list(),
                        mix_skurts =  list(), mix_fifths =  list(),
                        mix_sixths =  list(), marginal = list(),
                        lam  =  NULL, p_zip = 0, size = NULL,
                        prob = NULL, mu = NULL, p_zinb = 0, rho = NULL) {
  if (!is.null(Y_cat)) k_cat <- ncol(Y_cat) else k_cat <- 0
  if (!is.null(Y_cont)) k_cont <- ncol(Y_cont) else k_cont <- 0
  if (!is.null(Y_mix)) k_mix <- ncol(Y_mix) else k_mix <- 0
  if (!is.null(Y_pois)) k_pois <- ncol(Y_pois) else k_pois <- 0
  if (!is.null(Y_nb)) k_nb <- ncol(Y_nb) else k_nb <- 0
  if (k_pois > 0) {
    if (length(p_zip) < k_pois)
      p_zip <- c(rep(0, k_pois - length(p_zip)), p_zip)
  }
  if (k_nb > 0) {
    if (length(prob) > 0)
      mu <- size * (1 - prob)/prob
    if (length(p_zinb) < k_nb)
      p_zinb <- c(rep(0, k_nb - length(p_zinb)), p_zinb)
  }
  if (k_mix > 0) {
    if (class(mix_pis) != "list") {
      mix_pis <- list(mix_pis)
      mix_mus <- list(mix_mus)
      mix_sigmas <- list(mix_sigmas)
      mix_skews <- list(mix_skews)
      mix_skurts <- list(mix_skurts)
      if (length(mix_fifths) > 0) {
        mix_fifths <- list(mix_fifths)
        mix_sixths <- list(mix_sixths)
      }
    }
  }
  if (!is.null(rho)) {
    rho_calc <- cor(cbind(Y_cat, Y_cont, Y_comp, Y_pois, Y_nb))
    emax <- max(abs(rho_calc - rho))
  }
  result <- list()
  if (k_cat > 0) {
    ord_sum <- list()
    for (i in 1:k_cat) {
      ord_sum[[i]] <- as.data.frame(cbind(append(marginal[[i]], 1),
        cumsum(table(Y_cat[, i]))/n))
      colnames(ord_sum[[i]]) <- c("Target", "Simulated")
    }
    result <- append(result, list(ord_sum = ord_sum))
  }
  if ((k_cont + k_mix) > 0) {
    Yb <- cbind(Y_cont, Y_comp)
    cont_sum <- describe(Yb, type = 1)
    mom <- apply(Yb, 2, calc_moments)
    sim_fifths <- mom[5, ]
    sim_sixths <- mom[6, ]
    cont_sum <- as.data.frame(cbind(c(1:ncol(Yb)),
      cont_sum[, -c(1, 6, 7, 10, 13)], sim_fifths, sim_sixths))
    colnames(cont_sum) <- c("Distribution", "N", "Mean", "SD", "Median",
      "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
    means2 <- NULL
    vars2 <- NULL
    if (k_cont > 0) {
      means2 <- means[1:k_cont]
      vars2 <- vars[1:k_cont]
    }
    if (k_mix > 0) {
      means2 <- c(means2, unlist(mix_mus))
      vars2 <- c(vars2, (unlist(mix_sigmas))^2)
    }
    if ((length(fifths) == 0 & k_cont > 0) |
        (length(mix_fifths) == 0 & k_mix > 0)) {
      target_sum <- as.data.frame(cbind(c(1:ncol(Yb)), means2, vars2,
        c(skews, unlist(mix_skews)), c(skurts, unlist(mix_skurts))))
      colnames(target_sum) <- c("Distribution", "Mean", "SD", "Skew",
        "Skurtosis")
    } else {
      target_sum <- as.data.frame(cbind(c(1:ncol(Yb)), means2, vars2,
        c(skews, unlist(mix_skews)), c(skurts, unlist(mix_skurts)),
        c(fifths, unlist(mix_fifths)), c(sixths, unlist(mix_sixths))))
      colnames(target_sum) <- c("Distribution", "Mean", "SD", "Skew",
        "Skurtosis", "Fifth", "Sixth")
    }
    rownames(cont_sum) <- 1:ncol(Yb)
    rownames(target_sum) <- 1:ncol(Yb)
    result <- append(result, list(cont_sum = cont_sum,
      target_sum = target_sum))
    if (k_mix > 0) {
      target_mix <- NULL
      mix_sum <- describe(Y_mix, type = 1)
      mom <- apply(Y_mix, 2, calc_moments)
      sim_fifths <- mom[5, ]
      sim_sixths <- mom[6, ]
      mix_sum <- as.data.frame(cbind(c(1:k_mix),
        mix_sum[, -c(1, 6, 7, 10, 13)], sim_fifths, sim_sixths))
      colnames(mix_sum) <- c("Distribution", "N", "Mean", "SD", "Median",
        "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
      if (length(mix_fifths) == 0) {
        for (i in 1:k_mix) {
          target_mix <- rbind(target_mix, calc_mixmoments(mix_pis[[i]],
            mix_mus[[i]], mix_sigmas[[i]], mix_skews[[i]], mix_skurts[[i]]))
        }
        target_mix <- as.data.frame(cbind(c(1:k_mix), target_mix))
        colnames(target_mix) <- c("Distribution", "Mean", "SD", "Skew",
          "Skurtosis")
      } else {
        for (i in 1:k_mix) {
          target_mix <- rbind(target_mix, calc_mixmoments(mix_pis[[i]],
            mix_mus[[i]], mix_sigmas[[i]], mix_skews[[i]], mix_skurts[[i]],
            mix_fifths[[i]], mix_sixths[[i]]))
        }
        target_mix <- as.data.frame(cbind(c(1:k_mix), target_mix))
        colnames(target_mix) <- c("Distribution", "Mean", "SD", "Skew",
          "Skurtosis", "Fifth", "Sixth")
      }
      rownames(target_mix) <- 1:k_mix
      rownames(mix_sum) <- 1:k_mix
      result <- append(result, list(mix_sum = mix_sum,
        target_mix = target_mix))
      if (!is.null(rho)) {
        rho_mix <- cor(cbind(Y_cat, Y_cont, Y_mix, Y_pois, Y_nb))
        result <- append(result, list(rho_mix = rho_mix))
      }
    }
  }
  if (k_pois > 0) {
    pois_sum <- describe(Y_pois, type = 1)
    p_0 <- apply(Y_pois, 2, function(x) sum(x == 0)/n)
    pois_sum <- as.data.frame(cbind(pois_sum$vars, pois_sum$n,
      p_0, p_zip + (1 - p_zip) * exp(-lam), pois_sum$mean,
      (1 - p_zip) * lam, (pois_sum[, 4])^2,
      lam + (lam^2) * p_zip/(1 - p_zip), pois_sum$median,
      pois_sum$min, pois_sum$max, pois_sum$skew,
      pois_sum$kurtosis))
    colnames(pois_sum) <- c("Distribution", "N", "P0", "Exp_P0", "Mean",
      "Exp_Mean", "Var", "Exp_Var", "Median", "Min", "Max", "Skew",
      "Skurtosis")
    result <- append(result, list(pois_sum = pois_sum))
  }
  if (k_nb > 0) {
    nb_sum <- describe(Y_nb, type = 1)
    prob <- size/(mu + size)
    p_0 <- apply(Y_nb, 2, function(x) sum(x == 0)/n)
    nb_sum <- as.data.frame(cbind(nb_sum$vars, nb_sum$n, p_0,
      p_zinb + (1 - p_zinb) * (prob^size), prob, nb_sum$mean,
      (1 - p_zinb) * mu, (nb_sum[, 4])^2,
      (1 - p_zinb) * mu * (1 + mu * (p_zinb + 1/size)),
      nb_sum$median, nb_sum$min, nb_sum$max, nb_sum$skew,
      nb_sum$kurtosis))
    colnames(nb_sum) <- c("Distribution", "N", "P0", "Exp_P0",
      "Prob", "Mean", "Exp_Mean", "Var", "Exp_Var", "Median", "Min",
      "Max", "Skew", "Skurtosis")
    result <- append(result, list(nb_sum = nb_sum))
  }
  if (!is.null(rho))
    result <- append(result, list(rho_calc = rho_calc, maxerr = emax))
  result
}
