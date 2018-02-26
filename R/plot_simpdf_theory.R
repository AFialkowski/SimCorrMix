#' @title Plot Simulated Probability Density Function and Target PDF by Distribution Name or Function for Continuous or Count Variables
#'
#' @description This plots the PDF of simulated continuous or count (regular or zero-inflated, Poisson or Negative Binomial) data and
#'     overlays the target PDF (if \code{overlay} = TRUE), which is specified by distribution name (plus up to 4 parameters) or PDF
#'     function \code{fx} (plus support bounds).  If a continuous target distribution is provided (\code{cont_var = TRUE}), the simulated
#'     data \eqn{y} is scaled and then transformed (i.e. \eqn{y = sigma * scale(y) + mu}) so that it has the same mean (\eqn{mu}) and
#'     variance (\eqn{sigma^2}) as the target distribution.  The PDF's of continuous variables are shown as lines (using
#'     \code{\link[ggplot2]{geom_density}} and \code{ggplot2::geom_line}).  It works for valid or invalid power method PDF's.
#'     The PMF's of count variables are shown as vertical bar graphs (using \code{ggplot2::geom_col}).  The function returns a
#'     \code{\link[ggplot2]{ggplot2-package}} object so the user can save it or modify it as necessary.  The graph parameters
#'     (i.e. \code{title}, \code{sim_color}, \code{sim_lty}, \code{sim_size}, \code{target_color}, \code{target_lty}, \code{target_size},
#'     \code{legend.position}, \code{legend.justification}, \code{legend.text.size}, \code{title.text.size},
#'     \code{axis.text.size}, and \code{axis.title.size}) are inputs to the \code{\link[ggplot2]{ggplot2-package}} functions so information about
#'     valid inputs can be obtained from that package's documentation.
#' @param sim_y a vector of simulated data
#' @param title the title for the graph (default = "Simulated Probability Density Function")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value) on the x-axis
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value) on the x-axis
#' @param sim_color the line color for the simulated PDF (or column fill color in the case of
#'     \code{Dist} = "Poisson" or "Negative_Binomial")
#' @param sim_lty the line type for the simulated PDF (default = 1, solid line)
#' @param sim_size the line width for the simulated PDF
#' @param col_width width of column for simulated/target PMF of count variables (default = 0.5)
#' @param overlay if TRUE (default), the target distribution is also plotted given either a distribution name (and parameters)
#'     or PDF function fx (with bounds = ylower, yupper)
#' @param cont_var TRUE (default) for continuous variables, FALSE for count variables
#' @param target_color the line color for the target PDF (or column fill color in the case of
#'     \code{Dist} = "Poisson" or "Negative_Binomial")
#' @param target_lty the line type for the target PDF (default = 2, dashed line)
#' @param target_size the line width for the target PDF
#' @param Dist name of the distribution. The possible values are: "Benini", "Beta", "Beta-Normal", "Birnbaum-Saunders", "Chisq",
#'     "Exponential", "Exp-Geometric", "Exp-Logarithmic", "Exp-Poisson", "F", "Fisk", "Frechet", "Gamma", "Gaussian", "Gompertz",
#'     "Gumbel", "Kumaraswamy", "Laplace", "Lindley", "Logistic", \cr"Loggamma", "Lognormal", "Lomax", "Makeham", "Maxwell",
#'     "Nakagami", "Paralogistic", "Pareto", "Perks", "Rayleigh", "Rice", "Singh-Maddala", \cr"Skewnormal", "t", "Topp-Leone", "Triangular",
#'     "Uniform", "Weibull", "Poisson", and "Negative_Binomial".
#'     Please refer to the documentation for each package (either \code{\link[stats]{stats-package}}, \code{\link[VGAM]{VGAM-package}}, or
#'     \code{\link[triangle]{triangle}}) for information on appropriate parameter inputs.
#' @param params a vector of parameters (up to 4) for the desired distribution (keep NULL if \code{fx} supplied instead); for
#'     Poisson variables, must be lambda (mean) and the probability of a structural zero (use 0 for regular Poisson variables); for
#'     Negative Binomial variables, must be size, mean and the probability of a structural zero (use 0 for regular NB variables)
#' @param fx a PDF input as a function of x only, i.e. \code{fx = function(x) 0.5 * (x - 1)^2}; must return a scalar
#'     (keep NULL if \code{Dist} supplied instead)
#' @param lower the lower support bound for \code{fx}
#' @param upper the upper support bound for \code{fx}
#' @param legend.position the position of the legend
#' @param legend.justification the justification of the legend
#' @param legend.text.size the size of the legend labels
#' @param title.text.size the size of the plot title
#' @param axis.text.size the size of the axes text (tick labels)
#' @param axis.title.size the size of the axes titles
#' @import ggplot2
#' @importFrom VGAM dbenini rbenini dbetanorm rbetanorm dbisa rbisa ddagum rdagum dexpgeom rexpgeom dexplog rexplog
#'     dexppois rexppois dfisk rfisk dfrechet rfrechet dgompertz rgompertz dgumbel rgumbel dkumar rkumar dlaplace rlaplace dlind rlind
#'     dlgamma rlgamma dlomax rlomax dmakeham rmakeham dmaxwell rmaxwell dnaka rnaka dparalogistic
#'     rparalogistic dpareto rpareto dperks rperks dgenray rgenray drice rrice dsinmad rsinmad dskewnorm rskewnorm
#'     dtopple rtopple dzipois rzipois dzinegbin rzinegbin
#' @importFrom triangle dtriangle rtriangle
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @export
#' @keywords plot
#' @seealso \code{\link[SimMultiCorrData]{calc_theory}}, \code{\link[ggplot2]{ggplot}}
#' @return A \code{\link[ggplot2]{ggplot2-package}} object.
#' @references Please see the references for \code{\link[SimCorrMix]{plot_simtheory}}.
#'
#' @examples
#' # Using normal mixture variable from contmixvar1 example
#' Nmix <- contmixvar1(n = 1000, "Polynomial", means = 0, vars = 1,
#'   mix_pis = c(0.4, 0.6), mix_mus = c(-2, 2), mix_sigmas = c(1, 1),
#'   mix_skews = c(0, 0), mix_skurts = c(0, 0), mix_fifths = c(0, 0),
#'   mix_sixths = c(0, 0))
#' plot_simpdf_theory(Nmix$Y_mix[, 1],
#'   title = "Mixture of Normal Distributions",
#'   fx = function(x) 0.4 * dnorm(x, -2, 1) + 0.6 * dnorm(x, 2, 1),
#'   lower = -5, upper = 5)
#' \dontrun{
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
#' plot_simpdf_theory(Bmix$Y_mix[, 1], title = "Mixture of Beta Distributions",
#'   fx = function(x) mix_pis[1] * dbeta(x, 6, 3) + mix_pis[2] *
#'     dbeta(x, 4, 1.5) + mix_pis[3] * dbeta(x, 10, 20), lower = 0, upper = 1)
#' }
#'
plot_simpdf_theory <-
  function(sim_y, title = "Simulated Probability Density Function",
           ylower = NULL, yupper = NULL, sim_color = "dark blue", sim_lty = 1,
           sim_size = 1, col_width = 0.5, overlay = TRUE, cont_var = TRUE,
           target_color = "dark green", target_lty = 2, target_size = 1,
           Dist = c("Benini", "Beta", "Beta-Normal", "Birnbaum-Saunders",
           "Chisq", "Dagum", "Exponential", "Exp-Geometric", "Exp-Logarithmic",
           "Exp-Poisson", "F", "Fisk", "Frechet", "Gamma", "Gaussian",
           "Gompertz", "Gumbel", "Kumaraswamy", "Laplace", "Lindley",
           "Logistic", "Loggamma", "Lognormal", "Lomax", "Makeham", "Maxwell",
           "Nakagami", "Paralogistic", "Pareto", "Perks", "Rayleigh", "Rice",
           "Singh-Maddala", "Skewnormal", "t", "Topp-Leone", "Triangular",
           "Uniform", "Weibull", "Poisson", "Negative_Binomial"),
           params = NULL, fx = NULL, lower = NULL, upper = NULL,
           legend.position = c(0.975, 0.9), legend.justification = c(1, 1),
           legend.text.size = 10, title.text.size = 15, axis.text.size = 10,
           axis.title.size = 13) {
  if (overlay == FALSE) {
    if (is.null(ylower) & is.null(yupper)) {
      ylower <- min(sim_y)
      yupper <- max(sim_y)
    }
    if (cont_var == FALSE) {
      data <- as.data.frame(table(as.factor(sim_y))/length(sim_y))
      colnames(data) <- c("x", "y")
      data$type <- as.factor(rep("sim", nrow(data)))
      plot1 <- ggplot() + theme_bw() + ggtitle(title) +
        geom_col(data = data[data$type == "sim", ],
          width = col_width, aes_(x = ~x, y = ~y, fill = ~type),
          na.rm = TRUE) +
        xlab("y") + ylab("Probability") +
        theme(plot.title = element_text(size = title.text.size, face = "bold",
                                        hjust = 0.5),
              axis.text.x = element_text(size = axis.text.size),
              axis.title.x = element_text(size = axis.title.size),
              axis.text.y = element_text(size = axis.text.size),
              axis.title.y = element_text(size = axis.title.size),
              legend.text = element_text(size = legend.text.size),
              legend.position = legend.position,
              legend.justification = legend.justification) +
        scale_fill_manual(name = "", values = sim_color,
                          labels = c("Simulated Variable"))
    } else {
      data <- data.frame(x = 1:length(sim_y), y = sim_y,
                         type = as.factor(rep("sim", length(sim_y))))
      plot1 <- ggplot() + theme_bw() + ggtitle(title) +
        geom_density(data = data, aes_(x = ~y, colour = "Density",
          lty = ~type, size = ~type), na.rm = TRUE) +
        scale_x_continuous(name = "y", limits = c(ylower, yupper)) +
        scale_y_continuous(name = "Probability") +
        theme(plot.title = element_text(size = title.text.size, face = "bold",
                                        hjust = 0.5),
              axis.text.x = element_text(size = axis.text.size),
              axis.title.x = element_text(size = axis.title.size),
              axis.text.y = element_text(size = axis.text.size),
              axis.title.y = element_text(size = axis.title.size),
              legend.text = element_text(size = legend.text.size),
              legend.position = legend.position,
              legend.justification = legend.justification) +
        scale_linetype_manual(name = "", values = c(sim_lty),
                              labels = "Simulated Variable") +
        scale_colour_manual(name = "", values = c(sim_color),
                            labels = "Simulated Variable") +
        scale_size_manual(name = "", values = c(sim_size),
                            labels = "Simulated Variable")
    }
    return(plot1)
  }
  if (overlay == TRUE) {
    if (cont_var == TRUE) {
      if (!is.null(fx)) {
        theory <- calc_theory(fx = fx, lower = lower, upper = upper)
      }
      if (is.null(fx)) {
        theory <- calc_theory(Dist = Dist, params = params)
      }
      mu <- theory[1]
      sigma <- theory[2]
      sim_y <- sigma * scale(sim_y) + mu
    }
    if (is.null(ylower) & is.null(yupper)) {
      ylower <- min(sim_y)
      yupper <- max(sim_y)
    }
    x <- sim_y
    y_fx <- numeric(length(x))
    if (!is.null(fx)) {
      for (j in 1:length(x)) {
        y_fx[j] <- fx(x[j])
      }
    }
    if (is.null(fx)) {
      D <-
      data.frame(Dist = c("Benini", "Beta", "Beta-Normal", "Birnbaum-Saunders",
                          "Chisq", "Dagum", "Exponential", "Exp-Geometric",
                          "Exp-Logarithmic", "Exp-Poisson", "F", "Fisk",
                          "Frechet", "Gamma", "Gaussian", "Gompertz", "Gumbel",
                          "Kumaraswamy", "Laplace", "Lindley", "Logistic",
                          "Loggamma", "Lognormal", "Lomax",
                          "Makeham", "Maxwell", "Nakagami", "Paralogistic",
                          "Pareto", "Perks", "Rayleigh", "Rice",
                          "Singh-Maddala", "Skewnormal", "t", "Topp-Leone",
                          "Triangular", "Uniform", "Weibull", "Poisson",
                          "Negative_Binomial"),
                 pdf = c("dbenini", "dbeta", "dbetanorm", "dbisa", "dchisq",
                         "ddagum", "dexp", "dexpgeom", "dexplog", "dexppois",
                         "df", "dfisk", "dfrechet", "dgamma", "dnorm",
                         "dgompertz", "dgumbel", "dkumar", "dlaplace",
                         "dlind", "dlogis", "dlgamma", "dlnorm",
                         "dlomax", "dmakeham", "dmaxwell", "dnaka",
                         "dparalogistic", "dpareto", "dperks", "dgenray",
                         "drice", "dsinmad", "dskewnorm", "dt", "dtopple",
                         "dtriangle", "dunif", "dweibull", "dzipois",
                         "dzinegbin"),
                 fx = c("rbenini", "rbeta", "rbetanorm", "rbisa", "rchisq",
                        "rdagum", "rexp", "rexpgeom", "rexplog", "rexppois",
                        "rf", "rfisk", "rfrechet", "rgamma", "rnorm",
                        "rgompertz", "rgumbel", "rkumar", "rlaplace",
                        "rlind", "rlogis", "rlgamma", "rlnorm",
                        "rlomax", "rmakeham", "rmaxwell", "rnaka",
                        "rparalogistic", "rpareto", "rperks", "rgenray",
                        "rrice", "rsinmad", "rskewnorm", "rt", "rtopple",
                        "rtriangle", "runif", "rweibull", "rzipois",
                        "rzinegbin"),
                 Lower = as.numeric(c(params[1], 0, -Inf, rep(0, 9),
                                      params[1], 0, -Inf, 0, -Inf, 0, -Inf,
                                      0, -Inf, -Inf, rep(0, 6),
                                      params[1], rep(0, 4), -Inf, -Inf, 0,
                                      params[1], params[1], 0, 0, 0)),
                 Upper = as.numeric(c(Inf, 1, rep(Inf, 15), 1, rep(Inf, 17),
                                      1, params[2], params[2], Inf, Inf, Inf)))
      if (Dist == "Negative_Binomial") {
        y_fx <- dzinegbin(unique(sort(sim_y)), size = params[1],
                          munb = params[2], pstr0 = params[3])
      }
      if (Dist == "Poisson") {
        y_fx <- dzipois(unique(sort(sim_y)), params[1], params[2])
      }
      if (cont_var == TRUE) {
        i <- match(Dist, D$Dist)
        p <- as.character(D$pdf[i])
        if (length(params) == 1) y_fx <- get(p)(x, params[1])
        if (length(params) == 2) y_fx <- get(p)(x, params[1], params[2])
        if (length(params) == 3) y_fx <- get(p)(x, params[1], params[2],
                                                params[3])
        if (length(params) == 4) y_fx <- get(p)(x, params[1], params[2],
                                                params[3], params[4])
      }
    }
    if (cont_var == FALSE) {
      data <- as.data.frame(table(as.factor(sim_y))/length(sim_y))
      colnames(data) <- c("x", "y")
      data$type <- as.factor(rep("sim", nrow(data)))
      data2 <- data.frame(x = data$x, y = y_fx,
                          type = as.factor(rep("theory", length(y_fx))))
      data2 <- data.frame(rbind(data, data2))
      plot1 <- ggplot() + theme_bw() + ggtitle(title) +
        geom_col(data = data2[data2$type == "theory", ],
          width = col_width, aes_(x = ~x, y = ~y, fill = ~type),
          na.rm = TRUE) +
        geom_col(data = data2[data2$type == "sim", ],
          width = col_width, aes_(x = ~x, y = ~y, fill = ~type),
          na.rm = TRUE) +
        xlab("y") + ylab("Probability") +
        theme(plot.title = element_text(size = title.text.size, face = "bold",
                                        hjust = 0.5),
              axis.text.x = element_text(size = axis.text.size),
              axis.title.x = element_text(size = axis.title.size),
              axis.text.y = element_text(size = axis.text.size),
              axis.title.y = element_text(size = axis.title.size),
              legend.text = element_text(size = legend.text.size),
              legend.position = legend.position,
              legend.justification = legend.justification) +
        scale_fill_manual(name = "", values = c(sim_color, target_color),
                          labels = c("Simulated Variable", "Target"))
    } else {
      data <- data.frame(x = 1:length(sim_y), y = sim_y,
                         type = as.factor(rep("sim", length(sim_y))))
      data2 <- data.frame(x = x, y = y_fx,
                          type = as.factor(rep("theory", length(y_fx))))
      data2 <- data.frame(rbind(data, data2))
      plot1 <- ggplot() + theme_bw() + ggtitle(title) +
        geom_density(data = data2[data2$type == "sim", ],
          aes_(x = ~y, colour = ~type, lty = ~type, size = ~type),
          na.rm = TRUE) +
        geom_line(data = data2[data2$type == "theory", ],
          aes_(x = ~x, y = ~y, colour = ~type, lty = ~type, size = ~type),
          na.rm = TRUE) +
        scale_x_continuous(name = "y", limits = c(ylower, yupper)) +
        scale_y_continuous(name = "Probability") +
        theme(plot.title = element_text(size = title.text.size, face = "bold",
                                        hjust = 0.5),
              axis.text.x = element_text(size = axis.text.size),
              axis.title.x = element_text(size = axis.title.size),
              axis.text.y = element_text(size = axis.text.size),
              axis.title.y = element_text(size = axis.title.size),
              legend.text = element_text(size = legend.text.size),
              legend.position = legend.position,
              legend.justification = legend.justification) +
        scale_linetype_manual(name = "", values = c(sim_lty, target_lty),
                              labels = c("Simulated Variable", "Target")) +
        scale_colour_manual(name = "", values = c(sim_color, target_color),
                            labels = c("Simulated Variable", "Target")) +
        scale_size_manual(name = "", values = c(sim_size, target_size),
                            labels = c("Simulated Variable", "Target"))
    }
    return(plot1)
  }
}
