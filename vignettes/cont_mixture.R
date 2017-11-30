## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, fig.width = 6, fig.height = 4)

## ---- include=FALSE------------------------------------------------------
library("SimCorrMix")
library("printr")

## ------------------------------------------------------------------------
options(scipen = 999)
n <- 10000
mix_pis <- c(0.4, 0.6)
mix_mus <- c(-2, 2)
mix_sigmas <- c(1, 1)
mix_skews <- rep(0, 2)
mix_skurts <- rep(0, 2)
mix_fifths <- rep(0, 2)
mix_sixths <- rep(0, 2)
Nstcum <- calc_mixmoments(mix_pis, mix_mus, mix_sigmas, mix_skews, 
  mix_skurts, mix_fifths, mix_sixths)

## ------------------------------------------------------------------------
validpar(k_mix = 1, method = "Polynomial", means = Nstcum[1], 
  vars = Nstcum[2]^2, mix_pis = mix_pis, mix_mus = mix_mus, 
  mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts, 
  mix_fifths = mix_fifths, mix_sixths = mix_sixths)
Nmix2 <- contmixvar1(n, "Polynomial", Nstcum[1], Nstcum[2]^2, mix_pis, mix_mus, 
  mix_sigmas, mix_skews, mix_skurts, mix_fifths, mix_sixths)

## ------------------------------------------------------------------------
knitr::kable(Nmix2$target_mix, digits = 5, row.names = FALSE, 
  caption = "Summary of Target Distribution")
knitr::kable(Nmix2$mix_sum, digits = 5, row.names = FALSE, 
  caption = "Summary of Simulated Distribution")

## ------------------------------------------------------------------------
Nmix2$constants
Nmix2$valid.pdf

## ------------------------------------------------------------------------
fx <- function(x) 0.4 * dnorm(x, -2, 1) + 0.6 * dnorm(x, 2, 1)
cfx <- function(x, alpha, FUN = fx) {
  integrate(function(x, FUN = fx) FUN(x), -Inf, x, subdivisions = 1000, 
    stop.on.error = FALSE)$value - (1 - alpha)
}
y_star <- uniroot(cfx, c(3.3, 3.4), tol = 0.001, alpha = 0.05)$root
y_star

## ------------------------------------------------------------------------
sim_cdf_prob(sim_y = Nmix2$Y_mix[, 1], delta = y_star)$cumulative_prob

## ------------------------------------------------------------------------
plot_simpdf_theory(sim_y = Nmix2$Y_mix[, 1], ylower = -10, yupper = 10, 
  title = "PDF of Mixture of Normal Distributions", fx = fx, lower = -Inf, 
  upper = Inf)

## ------------------------------------------------------------------------
plot_sim_cdf(sim_y = Nmix2$Y_mix[, 1], calc_cprob = TRUE, delta = y_star)

## ------------------------------------------------------------------------
Stcum1 <- calc_theory("Beta", c(6, 3))
Stcum2 <- calc_theory("Beta", c(4, 1.5))
Stcum3 <- calc_theory("Beta", c(10, 20))
mix_pis <- c(0.5, 0.2, 0.3)
mix_mus <- c(Stcum1[1], Stcum2[1], Stcum3[1])
mix_sigmas <- c(Stcum1[2], Stcum2[2], Stcum3[2])
mix_skews <- c(Stcum1[3], Stcum2[3], Stcum3[3])
mix_skurts <- c(Stcum1[4], Stcum2[4], Stcum3[4])
mix_fifths <- c(Stcum1[5], Stcum2[5], Stcum3[5])
mix_sixths <- c(Stcum1[6], Stcum2[6], Stcum3[6])
mix_Six <- list(seq(0.01, 10, 0.01), c(0.01, 0.02, 0.03), seq(0.01, 10, 0.01))
Bstcum <- calc_mixmoments(mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, 
  mix_fifths, mix_sixths)
validpar(k_mix = 1, method = "Polynomial", means = Bstcum[1], 
  vars = Bstcum[2]^2, mix_pis = mix_pis, mix_mus = mix_mus, 
  mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts, 
  mix_fifths = mix_fifths, mix_sixths = mix_sixths, mix_Six = mix_Six)
Bmix3 <- contmixvar1(n, "Polynomial", Bstcum[1], Bstcum[2]^2, mix_pis, mix_mus, 
  mix_sigmas, mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six)
knitr::kable(Bmix3$target_mix, digits = 5, row.names = FALSE, 
  caption = "Summary of Target Distribution")
knitr::kable(Bmix3$mix_sum, digits = 5, row.names = FALSE, 
  caption = "Summary of Simulated Distribution")

## ------------------------------------------------------------------------
Bmix3$valid.pdf
Bmix3$sixth_correction
Bmix3$constants
Bplot <- plot_simpdf_theory(Bmix3$Y_mix[, 1], 
  title = "PDF of Mixture of Beta Distributions",
  fx = function(x) mix_pis[1] * dbeta(x, 6, 3) + 
    mix_pis[2] * dbeta(x, 4, 1.5) + mix_pis[3] * dbeta(x, 10, 20), lower = 0, 
  upper = 1)
Bplot

## ------------------------------------------------------------------------
Stcum1 <- calc_theory("Chisq", 2)
Stcum2 <- calc_theory("Chisq", 32)
mix_pis <- c(0.7, 0.3)
mix_mus <- c(Stcum1[1], Stcum2[1])
mix_sigmas <- c(Stcum1[2], Stcum2[2])
mix_skews <- c(Stcum1[3], Stcum2[3])
mix_skurts <- c(Stcum1[4], Stcum2[4])
mix_fifths <- c(Stcum1[5], Stcum2[5])
mix_sixths <- c(Stcum1[6], Stcum2[6])
Cstcum <- calc_mixmoments(mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, 
  mix_fifths, mix_sixths)
validpar(k_mix = 1, method = "Polynomial", means = Cstcum[1], 
  vars = Cstcum[2]^2, mix_pis = mix_pis, mix_mus = mix_mus, 
  mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts, 
  mix_fifths = mix_fifths, mix_sixths = mix_sixths)
Cmix2 <- contmixvar1(n, "Polynomial", Cstcum[1], Cstcum[2]^2, mix_pis, mix_mus, 
  mix_sigmas, mix_skews, mix_skurts, mix_fifths, mix_sixths)
knitr::kable(Cmix2$target_mix, digits = 5, row.names = FALSE, 
  caption = "Summary of Target Distribution")
knitr::kable(Cmix2$mix_sum, digits = 5, row.names = FALSE, 
  caption = "Summary of Simulated Distribution")

## ------------------------------------------------------------------------
Cmix2$valid.pdf
Cmix2$constants
Cplot <- plot_simpdf_theory(Cmix2$Y_mix[, 1], 
  title = "PDF of Mixture of Chi-square Distributions",
  fx = function(x) mix_pis[1] * dchisq(x, 2) + mix_pis[2] * dchisq(x, 32), 
  lower = 0, upper = Inf)
Cplot

## ------------------------------------------------------------------------
mix_pis <- c(0.3, 0.2, 0.4, 0.1)
mix_mus <- c(-4, -2, 1, 5)
mix_sigmas <- c(1, 2, 3, 4)
mix_skews <- rep(0, 4)
mix_skurts <- rep(0, 4)
mix_fifths <- rep(0, 4)
mix_sixths <- rep(0, 4)
Nstcum <- calc_mixmoments(mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, 
  mix_fifths, mix_sixths)
validpar(k_mix = 1, method = "Polynomial", means = Nstcum[1], 
  vars = Nstcum[2]^2, mix_pis = mix_pis, mix_mus = mix_mus, 
  mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts, 
  mix_fifths = mix_fifths, mix_sixths = mix_sixths)
Nmix4 <- contmixvar1(n, "Polynomial", Nstcum[1], Nstcum[2]^2, mix_pis, mix_mus, 
  mix_sigmas, mix_skews, mix_skurts, mix_fifths, mix_sixths)
knitr::kable(Nmix4$target_mix, digits = 5, row.names = FALSE, 
  caption = "Summary of Target Distribution")
knitr::kable(Nmix4$mix_sum, digits = 5, row.names = FALSE, 
  caption = "Summary of Simulated Distribution")

## ------------------------------------------------------------------------
Nmix4$valid.pdf
Nmix4$constants
Nplot <- plot_simpdf_theory(Nmix4$Y_mix[, 1], 
  title = "PDF of Mixture of Normal Distributions",
  fx = function(x) mix_pis[1] * dnorm(x, mix_mus[1], mix_sigmas[1]) + 
    mix_pis[2] * dnorm(x, mix_mus[2], mix_sigmas[2]) + 
    mix_pis[3] * dnorm(x, mix_mus[3], mix_sigmas[3]) + 
    mix_pis[4] * dnorm(x, mix_mus[4], mix_sigmas[4]), 
  lower = -Inf, upper = Inf)
Nplot

