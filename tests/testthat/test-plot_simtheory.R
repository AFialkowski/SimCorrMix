skip_on_cran()
library("SimCorrMix")
context("Plot simulated data")

options(scipen = 999)
tol <- 1e-5

seed <- 276
n <- 10000

# Continuous variables
L <- calc_theory("Logistic", c(0, 1))
C <- calc_theory("Chisq", 4)
B <- calc_theory("Beta", c(4, 1.5))

# Non-mixture variable
skews <- L[3]
skurts <- L[4]
fifths <- L[5]
sixths <- L[6]
Six <- list(1.75)

# Mixture variable
mix_pis <- list(c(0.3, 0.2, 0.5))
mix_mus <- list(c(L[1], C[1], B[1]))
mix_sigmas <- list(c(L[2], C[2], B[2]))
mix_skews <- list(c(L[3], C[3], B[3]))
mix_skurts <- list(c(L[4], C[4], B[4]))
mix_fifths <- list(c(L[5], C[5], B[5]))
mix_sixths <- list(c(L[6], C[6], B[6]))
mix_Six <- list(1.75, NULL, 0.03)
Mstcum <- calc_mixmoments(mix_pis[[1]], mix_mus[[1]], mix_sigmas[[1]],
  mix_skews[[1]], mix_skurts[[1]], mix_fifths[[1]], mix_sixths[[1]])

means <- c(L[1], Mstcum[1])
vars <- c(L[2]^2, Mstcum[2]^2)

marginal <- list(0.3)
support <- list(c(0, 1))
lam <- c(0.5, 5)
p_zip <- c(0, 0.1)
size <- c(2, 5)
prob <- c(0.75, 0.5)
mu <- size * (1 - prob)/prob
p_zinb <- c(0, 0.2)

k_cat <- length(marginal)
k_cont <- length(Six)
k_mix <- length(mix_pis)
k_comp <- sum(unlist(lapply(mix_pis, length)))
k_pois <- length(lam)
k_nb <- length(size)
k_total <- k_cat + k_cont + k_comp + k_pois + k_nb

Rey <- matrix(0.35, k_total, k_total)
diag(Rey) <- 1
rownames(Rey) <- colnames(Rey) <- c("O1", "C1", "M1_1", "M1_2", "M1_3", "P1",
  "ZIP1", "NB1", "ZINB1")
Rey["M1_1", "M1_2"] <- Rey["M1_2", "M1_1"] <- Rey["M1_1", "M1_3"] <-
  Rey["M1_3", "M1_1"] <- Rey["M1_2", "M1_3"] <- Rey["M1_3", "M1_2"] <- 0

Sim1 <- corrvar(n, k_cat, k_cont, k_mix, k_pois, k_nb,
  "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
  mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
  mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
  mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
  support = support, lam = lam, p_zip = p_zip[2], size = size,
  prob = prob, p_zinb = p_zinb[2], rho = Rey, seed = seed, epsilon = 0.01)

test_that("works for continuous variable", {
  expect_is(plot_simtheory(Sim1$Y_comp[, 2],
    overlay = FALSE), "ggplot")
  expect_is(plot_simtheory(Sim1$Y_comp[, 2],
    overlay = TRUE, Dist = "Chisq", params = 4), "ggplot")
})

test_that("works for Poisson variable", {
  expect_is(plot_simtheory(Sim1$Y_pois[, 1], cont_var = FALSE,
    overlay = FALSE), "ggplot")
  expect_is(plot_simtheory(Sim1$Y_pois[, 1], cont_var = FALSE,
    Dist = "Poisson", params = c(lam[1], p_zip[1])), "ggplot")
  expect_is(plot_simtheory(Sim1$Y_pois[, 2], cont_var = FALSE,
    overlay = FALSE), "ggplot")
  expect_is(plot_simtheory(Sim1$Y_pois[, 2], cont_var = FALSE,
    Dist = "Poisson", params = c(lam[2], p_zip[2])), "ggplot")
})

test_that("works for NB variable", {
  expect_is(plot_simtheory(Sim1$Y_nb[, 1], cont_var = FALSE,
    overlay = FALSE), "ggplot")
  expect_is(plot_simtheory(Sim1$Y_nb[, 1], cont_var = FALSE,
    Dist = "Negative_Binomial", params = c(size[1], mu[1], p_zinb[1])),
    "ggplot")
  expect_is(plot_simtheory(Sim1$Y_nb[, 2], cont_var = FALSE,
    overlay = FALSE), "ggplot")
  expect_is(plot_simtheory(Sim1$Y_nb[, 2], cont_var = FALSE,
    Dist = "Negative_Binomial", params = c(size[2], mu[2], p_zinb[2])),
    "ggplot")
})
