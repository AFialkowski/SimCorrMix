#' @title Plot Simulated Data and Target Distribution Data by Name or Function for Continuous or Count Variables
#'
#' @description This plots simulated continuous or count (regular or zero-inflated, Poisson or Negative Binomial) data and overlays data
#'     (if \code{overlay} = TRUE) generated from the target distribution.  The target is specified by name (plus up to 4 parameters) or
#'     PDF function \code{fx} (plus support bounds).  Due to the integration involved in finding the CDF from the PDF supplied by
#'     \code{fx}, only continuous \code{fx} may be supplied.  Both are plotted as histograms (using \code{\link[ggplot2]{geom_histogram}}).
#'     If a continuous target distribution is specified (\code{cont_var = TRUE}), the simulated data \eqn{y} is
#'     scaled and then transformed (i.e. \eqn{y = sigma * scale(y) + mu}) so that it has the same mean (\eqn{mu}) and variance
#'     (\eqn{sigma^2}) as the target distribution.  It works for valid or invalid power method PDF's.  It returns a
#'     \code{\link[ggplot2]{ggplot2-package}} object so the user can save it or
#'     modify it as necessary.  The graph parameters (i.e. \code{title}, \code{sim_color}, \code{target_color},
#'     \code{legend.position}, \code{legend.justification}, \code{legend.text.size}, \code{title.text.size},
#'     \code{axis.text.size}, and \code{axis.title.size}) are inputs to the \code{\link[ggplot2]{ggplot2-package}} functions so information about
#'     valid inputs can be obtained from that package's documentation.
#' @param sim_y a vector of simulated data
#' @param title the title for the graph (default = "Simulated Data Values")
#' @param ylower the lower y value to use in the plot (default = NULL, uses minimum simulated y value) on the y-axis
#' @param yupper the upper y value (default = NULL, uses maximum simulated y value) on the y-axis
#' @param sim_color the histogram fill color for the simulated variable (default = "dark blue")
#' @param overlay if TRUE (default), the target distribution is also plotted given either a distribution name (and parameters)
#'     or PDF function fx (with support bounds = lower, upper)
#' @param cont_var TRUE (default) for continuous variables, FALSE for count variables
#' @param target_color the histogram fill color for the target distribution (default = "dark green")
#' @param binwidth the width of bins to use when creating the histograms (default = NULL)
#' @param nbins the number of bins to use when creating the histograms (default = 100); overridden by \code{binwidth}
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
#' @param lower the lower support bound for a supplied \code{fx}, else keep NULL (note: if an error is thrown from \code{uniroot},
#'     try a slightly higher lower bound; i.e., 0.0001 instead of 0)
#' @param upper the upper support bound for a supplied \code{fx}, else keep NULL (note: if an error is thrown from \code{uniroot},
#'     try a lower upper bound; i.e., 100000 instead of Inf)
#' @param seed the seed value for random number generation (default = 1234)
#' @param sub the number of subdivisions to use in the integration to calculate the CDF from \code{fx}; if no result, try increasing
#'     sub (requires longer computation time; default = 1000)
#' @param legend.position the position of the legend
#' @param legend.justification the justification of the legend
#' @param legend.text.size the size of the legend labels
#' @param title.text.size the size of the plot title
#' @param axis.text.size the size of the axes text (tick labels)
#' @param axis.title.size the size of the axes titles
#' @import ggplot2
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var integrate
#' @import utils
#' @importFrom VGAM dbenini rbenini dbetanorm rbetanorm dbisa rbisa ddagum rdagum dexpgeom rexpgeom dexplog rexplog
#'     dexppois rexppois dfisk rfisk dfrechet rfrechet dgompertz rgompertz dgumbel rgumbel dkumar rkumar dlaplace rlaplace dlind rlind
#'     dlgamma rlgamma dlomax rlomax dmakeham rmakeham dmaxwell rmaxwell dnaka rnaka dparalogistic
#'     rparalogistic dpareto rpareto dperks rperks dgenray rgenray drice rrice dsinmad rsinmad dskewnorm rskewnorm
#'     dtopple rtopple dzipois rzipois dzinegbin rzinegbin
#' @importFrom triangle dtriangle rtriangle
#' @export
#' @keywords plot
#' @seealso \code{\link[SimMultiCorrData]{calc_theory}},
#'     \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_histogram}}
#' @return A \code{\link[ggplot2]{ggplot2-package}} object.
#' @references
#' Carnell R (2017). triangle: Provides the Standard Distribution Functions for the Triangle Distribution. R package version 0.11.
#'     \url{https://CRAN.R-project.org/package=triangle}.
#'
#' Fialkowski AC (2017). SimMultiCorrData: Simulation of Correlated Data with Multiple Variable Types. R package version 0.2.1.
#'     \url{https://CRAN.R-project.org/package=SimMultiCorrData}.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3):1-17. \cr \doi{10.18637/jss.v019.i03}.
#'
#' Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
#'
#' Yee TW (2017). VGAM: Vector Generalized Linear and Additive Models. \cr \url{https://CRAN.R-project.org/package=VGAM}.
#'
#' @examples
#' # Using normal mixture variable from contmixvar1 example
#' Nmix <- contmixvar1(n = 1000, "Polynomial", means = 0, vars = 1,
#'   mix_pis = c(0.4, 0.6), mix_mus = c(-2, 2), mix_sigmas = c(1, 1),
#'   mix_skews = c(0, 0), mix_skurts = c(0, 0), mix_fifths = c(0, 0),
#'   mix_sixths = c(0, 0))
#' plot_simtheory(Nmix$Y_mix[, 1], title = "Mixture of Normal Distributions",
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
#' plot_simtheory(Bmix$Y_mix[, 1], title = "Mixture of Beta Distributions",
#'   fx = function(x) mix_pis[1] * dbeta(x, 6, 3) + mix_pis[2] *
#'     dbeta(x, 4, 1.5) + mix_pis[3] * dbeta(x, 10, 20), lower = 0, upper = 1)
#' }
#'
plot_simtheory <- function(sim_y, title = "Simulated Data Values",
                           ylower = NULL, yupper = NULL,
                           sim_color = "dark blue", overlay = TRUE,
                           cont_var = TRUE, target_color = "dark green",
                           binwidth = NULL, nbins = 100,
                           Dist = c("Benini", "Beta", "Beta-Normal",
                                    "Birnbaum-Saunders", "Chisq", "Dagum",
                                    "Exponential", "Exp-Geometric",
                                    "Exp-Logarithmic", "Exp-Poisson", "F",
                                    "Fisk", "Frechet", "Gamma", "Gaussian",
                                    "Gompertz", "Gumbel", "Kumaraswamy",
                                    "Laplace", "Lindley", "Logistic",
                                    "Loggamma", "Lognormal", "Lomax",
                                    "Makeham", "Maxwell", "Nakagami",
                                    "Paralogistic", "Pareto", "Perks",
                                    "Rayleigh", "Rice", "Singh-Maddala",
                                    "Skewnormal", "t", "Topp-Leone",
                                    "Triangular", "Uniform", "Weibull",
                                    "Poisson", "Negative_Binomial"),
                           params = NULL, fx = NULL, lower = NULL,
                           upper = NULL, seed = 1234, sub = 1000,
                           legend.position = c(0.975, 0.9),
                           legend.justification = c(1, 1),
                           legend.text.size = 10, title.text.size = 15,
                           axis.text.size = 10, axis.title.size = 13) {
  if (overlay == FALSE) {
    if (is.null(ylower) & is.null(yupper)) {
      ylower <- min(sim_y)
      yupper <- max(sim_y)
    }
    data <- data.frame(x = 1:length(sim_y), y = sim_y,
                       type = as.factor(rep("sim", length(sim_y))))
    if (cont_var == FALSE)
      limits0 <- NULL else
        limits0 <- c(ylower, yupper)
    plot1 <- ggplot() + theme_bw() + ggtitle(title) +
      geom_histogram(data = data[data$type == "sim", ],
        aes_(~y, fill = ~type), binwidth = binwidth, bins = nbins,
        na.rm = TRUE) +
      scale_x_continuous(name = "y", limits = limits0) +
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
    set.seed(seed)
    n <- length(sim_y)
    if (!is.null(fx)) {
      uni <- runif(n)
      cfx <- function(x, u, FUN = fx) integrate(function(x, FUN = fx) FUN(x),
                                                lower, x, subdivisions = sub,
                                                stop.on.error = FALSE)$value -
        u
      y_fx <- rep(NA, n)
      for (i in 1:n) {
        y_fx[i] <- uniroot(cfx, c(lower, upper), extendInt = "yes",
                           tol = 0.0001, u = uni[i])$root
      }
      limits0 <- c(ylower, yupper)
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
      set.seed(seed)
      if (Dist == "Negative_Binomial") {
        y_fx <- rzinegbin(n, size = params[1], munb = params[2],
                          pstr0 = params[3])
      } else {
        i <- match(Dist, D$Dist)
        p <- as.character(D$fx[i])
        if (length(params) == 1) y_fx <- get(p)(n, params[1])
        if (length(params) == 2) y_fx <- get(p)(n, params[1], params[2])
        if (length(params) == 3) y_fx <- get(p)(n, params[1], params[2],
                                                params[3])
        if (length(params) == 4) y_fx <- get(p)(n, params[1], params[2],
                                                params[3], params[4])
      }
      if (Dist == "Poisson" | Dist == "Negative_Binomial")
        limits0 <- NULL else
          limits0 <- c(ylower, yupper)
    }
    data <- data.frame(x = 1:length(sim_y), y = sim_y,
                       type = as.factor(rep("sim", length(sim_y))))
    data2 <- data.frame(x = 1:length(y_fx), y = sort(y_fx),
                        type = as.factor(rep("theory", length(y_fx))))
    data2 <- data.frame(rbind(data, data2))
    plot1 <- ggplot() + theme_bw() + ggtitle(title) +
      geom_histogram(data = data2[data2$type == "theory", ],
        aes_(~y, fill = ~type), binwidth = binwidth, bins = nbins,
        na.rm = TRUE) +
      geom_histogram(data = data2[data2$type == "sim", ],
        aes_(~y, fill = ~type), binwidth = binwidth, bins = nbins,
        na.rm = TRUE) +
      scale_x_continuous(name = "y", limits = limits0) +
      theme(plot.title = element_text(size = title.text.size, face = "bold",
                                      hjust = 0.5),
            axis.text.x = element_text(size = axis.text.size),
            axis.title.x = element_text(size = axis.title.size),
            axis.text.y = element_text(size = axis.text.size),
            axis.title.y = element_text(size = axis.title.size),
            legend.text = element_text(size = legend.text.size),
            legend.position = legend.position,
            legend.justification = legend.justification) +
      scale_fill_manual(name = "", values = c("sim" = sim_color,
                                              "theory" = target_color),
                        labels = c("Simulated Variable", "Target Variable"))
    return(plot1)
  }
}
