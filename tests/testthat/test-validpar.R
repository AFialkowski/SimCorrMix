skip_on_cran()
library("SimCorrMix")
context("Check parameter inputs")

options(scipen = 999)
tol <- 1e-5

n <- 25
set.seed(1234)
cstart1 <- runif(n, min = -2, max = 2)
cstart2 <- runif(n, min = -1, max = 1)
cstart3 <- runif(n, min = -0.5, max = 0.5)
cstartF <- cbind(cstart1, cstart2, cstart3)
cstartF <- list(cstartF, cstartF, cstartF, cstartF)
cstartF2 <- cbind(cstart1, cstart2)
cstartF2 <- list(cstartF2, cstartF2, cstartF2, cstartF2)

set.seed(1234)
cstart1 <- runif(n, min = -2, max = 2)
cstart2 <- runif(n, min = -1, max = 1)
cstart3 <- runif(n, min = -1, max = 1)
cstart4 <- runif(n, min = -0.025, max = 0.025)
cstart5 <- runif(n, min = -0.025, max = 0.025)
cstartP <- cbind(cstart1, cstart2, cstart3, cstart4, cstart5)
cstartP <- list(cstartP, cstartP, cstartP, cstartP)
cstartP2 <- cbind(cstart1, cstart2, cstart3)
cstartP2 <- list(cstartP2, cstartP2, cstartP2, cstartP2)

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

test_that("works for ordinal variables", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = list(), lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey))
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = list(0.5, 0.2), lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey))
})

test_that("works to test length of method", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    method = c("Fleishman", "Polynomial"), means, vars, skews, skurts,
    mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey))
})

test_that("works to test length of means and vars", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means[1], vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey))
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars[2], skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey))
})

test_that("works to test length of skews", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews = c(0, 0), skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey))
})

test_that("works to test length of mix_skews", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = list(),
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey))
})

test_that("works to test length of Six", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths,
    Six = list(1.75, 1.75),
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six,
    marginal = marginal, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey))
})

test_that("works to test length of mix_Six", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six[1:2],
    marginal = marginal, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey))
})

test_that("works to test sum of mix_pis", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = list(c(0.3, 0.3, 0.5)), mix_mus = mix_mus,
    mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts,
    mix_fifths = mix_fifths, mix_sixths = mix_sixths, mix_Six = mix_Six,
    marginal = marginal, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey))
})

test_that("works to test length of cstart", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus,
    mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts,
    mix_fifths = mix_fifths, mix_sixths = mix_sixths, mix_Six = mix_Six,
    marginal = marginal, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey, cstart = cstartP[1:3]))
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey, cstart = cstartF[1:3]))
})

test_that("works to test dim of cstart", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus,
    mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts,
    mix_fifths = mix_fifths, mix_sixths = mix_sixths, mix_Six = mix_Six,
    marginal = marginal, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey, cstart = cstartP2))
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey, cstart = cstartF2))
})

test_that("works to test lam", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus,
    mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts,
    mix_fifths = mix_fifths, mix_sixths = mix_sixths, mix_Six = mix_Six,
    marginal = marginal, lam = lam[1], p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey))
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = c(0.5, -5),
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey))
})

test_that("works to test size, prob, and mu", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus,
    mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts,
    mix_fifths = mix_fifths, mix_sixths = mix_sixths, mix_Six = mix_Six,
    marginal = marginal, lam = lam, p_zip = p_zip[2], size = size[1],
    prob = prob, p_zinb = p_zinb[2], rho = Rey))
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob[1], p_zinb = p_zinb[2],
    rho = Rey))
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, mu = mu[1], p_zinb = p_zinb[2],
    rho = Rey))
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, mu = mu,
    p_zinb = p_zinb[2], rho = Rey))
})

test_that("works to test Sigma and rho", {
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus,
    mix_sigmas = mix_sigmas, mix_skews = mix_skews, mix_skurts = mix_skurts,
    mix_fifths = mix_fifths, mix_sixths = mix_sixths, mix_Six = mix_Six,
    marginal = marginal, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], Sigma = Rey[-1, 1], rho = Rey))
  expect_error(validpar(k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, lam = lam,
    p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey[-1, -1]))
})
