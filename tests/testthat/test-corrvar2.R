skip_on_cran()
library("SimCorrMix")
context("Simulate using correlation method 2")

options(scipen = 999)
tol <- 1e-5

n <- 25
set.seed(1234)
cstart1 <- runif(n, min = -2, max = 2)
cstart2 <- runif(n, min = -1, max = 1)
cstart3 <- runif(n, min = -0.5, max = 0.5)
cstartF <- cbind(cstart1, cstart2, cstart3)
cstartF <- list(cstartF, cstartF, cstartF, cstartF)

set.seed(1234)
cstart1 <- runif(n, min = -2, max = 2)
cstart2 <- runif(n, min = -1, max = 1)
cstart3 <- runif(n, min = -1, max = 1)
cstart4 <- runif(n, min = -0.025, max = 0.025)
cstart5 <- runif(n, min = -0.025, max = 0.025)
cstartP <- cbind(cstart1, cstart2, cstart3, cstart4, cstart5)
cstartP <- list(cstartP, cstartP, cstartP, cstartP)

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
pois_eps <- nb_eps <- 0.0001

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

test_that("works for Fleishman method: 1 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey, seed = seed, epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey, seed = seed, epsilon = 0.01,
    cstart = cstartF)$constants[1, "c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey, seed = seed, epsilon = 0.01,
    errorloop = TRUE)$constants[1, "c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 0 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat = 0, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey[-1, -1], seed = seed,
    epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat = 0, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, lam = lam, p_zip = p_zip[2], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-1, -1], seed = seed,
    epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat = 0, k_cont, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, lam = lam, p_zip = p_zip[2], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-1, -1], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 ordinal, 0 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Fleishman", means[2], vars[2], mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey[-2, -2], seed = seed, epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Fleishman", means[2], vars[2], mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey[-2, -2], seed = seed, epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Fleishman", means[2], vars[2], mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey[-2, -2], seed = seed, epsilon = 0.01,
    errorloop = TRUE)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 ordinal, 1 continuous, 0 mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix = 0, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, marginal = marginal,
    lam = lam, p_zip = p_zip, size = size, prob = prob,
    p_zinb = p_zinb, rho = Rey[-c(3:5), -c(3:5)], seed = seed,
    epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix = 0, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, marginal = marginal,
    lam = lam, p_zip = p_zip, size = size, mu = mu,
    p_zinb = p_zinb, rho = Rey[-c(3:5), -c(3:5)], seed = seed,
    epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix = 0, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts, marginal = marginal,
    lam = lam, p_zip = p_zip, size = size, prob = prob,
    p_zinb = p_zinb, rho = Rey[-c(3:5), -c(3:5)], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 ordinal, 0 continuous, 0 mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix = 0, k_pois, k_nb,
    marginal = marginal, support = support, lam = lam, p_zip = p_zip,
    size = size, prob = prob, p_zinb = p_zinb, rho = Rey[-c(2:5), -c(2:5)],
    seed = seed, epsilon = 0.001)$Sigma[1, 1], 1, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix = 0, k_pois, k_nb,
    marginal = marginal, support = support, lam = lam, p_zip = p_zip,
    size = size, mu = mu, p_zinb = p_zinb, rho = Rey[-c(2:5), -c(2:5)],
    seed = seed, epsilon = 0.001)$Sigma[1, 1], 1, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix = 0, k_pois, k_nb,
    marginal = marginal, support = support, lam = lam, p_zip = p_zip,
    size = size, mu = mu, p_zinb = p_zinb, rho = Rey[-c(2:5), -c(2:5)],
    seed = seed, epsilon = 0.001, errorloop = TRUE)$Sigma[1, 1], 1,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 ordinal, 1 continuous, 1 mixture,
          0 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam[2], p_zip = p_zip[2], size = size, prob = prob,
    p_zinb = p_zinb[2], rho = Rey[-6, -6], seed = seed,
    epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam[2], p_zip = p_zip[2], size = size, mu = mu,
    p_zinb = p_zinb[2], rho = Rey[-6, -6], seed = seed,
    epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam[2], p_zip = p_zip[2], size = size, mu = mu,
    p_zinb = p_zinb[2], rho = Rey[-6, -6], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 0 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam[1], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey[-7, -7], seed = seed, epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam[1], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey[-7, -7], seed = seed, epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam[1], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey[-7, -7], seed = seed, epsilon = 0.01,
    errorloop = TRUE)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 0 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size[2], prob = prob[2],
    p_zinb = p_zinb[2], rho = Rey[-8, -8], seed = seed,
    epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size[2], mu = mu[2],
    p_zinb = p_zinb[2], rho = Rey[-8, -8], seed = seed,
    epsilon = 0.01)$constants[1, "c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size[2], prob = prob[2],
    p_zinb = p_zinb[2], rho = Rey[-8, -8], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method: 1 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 1 NB, 0 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size[1], prob = prob[1],
    rho = Rey[-9, -9], seed = seed, epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size[1], mu = mu[1],
    rho = Rey[-9, -9], seed = seed, epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Fleishman", means, vars, skews, skurts, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size[1], prob = prob[1],
    rho = Rey[-9, -9], seed = seed, epsilon = 0.01,
    errorloop = TRUE)$constants[1, "c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

##################
## Polynomial
##################

test_that("works for Polynomial method: 1 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey, seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey, seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey, seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 0 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat = 0, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six,
    lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey[-1, -1], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat = 0, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six,
    lam = lam, p_zip = p_zip[2], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-1, -1], seed = seed,
    epsilon = 0.01, cstart = cstartP)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat = 0, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six,
    lam = lam, p_zip = p_zip[2], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-1, -1], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 0 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Polynomial", means[2], vars[2],
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip, size = size,
    prob = prob, p_zinb = p_zinb, rho = Rey[-2, -2], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Polynomial", means[2], vars[2],
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip, size = size,
    mu = mu, p_zinb = p_zinb, rho = Rey[-2, -2], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Polynomial", means[2], vars[2],
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip, size = size,
    prob = prob, p_zinb = p_zinb, rho = Rey[-2, -2], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 1 continuous, 0 mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix = 0, k_pois, k_nb,
    "Polynomial", means[1], vars[1], skews, skurts, fifths, sixths, Six,
    marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey[-c(3:5), -c(3:5)], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix = 0, k_pois, k_nb,
    "Polynomial", means[1], vars[1], skews, skurts, fifths, sixths, Six,
    marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-c(3:5), -c(3:5)], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix = 0, k_pois, k_nb,
    "Polynomial", means[1], vars[1], skews, skurts, fifths, sixths, Six,
    marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey[-c(3:5), -c(3:5)], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 1 continuous, 1 mixture,
          0 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam[2], p_zip = p_zip[2], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey[-6, -6], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam[2], p_zip = p_zip[2], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-6, -6], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam[2], p_zip = p_zip[2], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-6, -6], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 0 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam[1], size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey[-7, -7], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam[1], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-7, -7], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 1, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam[1], size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-7, -7], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 0 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size[2],
    prob = prob[2], p_zinb = p_zinb[2], rho = Rey[-8, -8], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size[2],
    mu = mu[2], p_zinb = p_zinb[2], rho = Rey[-8, -8], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size[2],
    prob = prob[2], p_zinb = p_zinb[2], rho = Rey[-8, -8], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 1 NB, 0 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size[1],
    prob = prob[1], rho = Rey[-9, -9], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size[1],
    mu = mu[1], rho = Rey[-9, -9], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 1,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2], size = size[1],
    mu = mu[1], rho = Rey[-9, -9], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 1 continuous, 1 mixture,
          0 Poisson, 0 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 0, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, size = size,
    prob = prob, p_zinb = p_zinb[2], rho = Rey[-c(6:7), -c(6:7)], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 0, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-c(6:7), -c(6:7)], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois = 0, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, size = size,
    mu = mu, p_zinb = p_zinb[2], rho = Rey[-c(6:7), -c(6:7)], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 1 continuous, 1 mixture,
          1 Poisson, 1 ZIP, 0 NB, 0 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 0,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2],
    rho = Rey[-c(8:9), -c(8:9)], seed = seed,
    epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb = 0,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    mix_pis = mix_pis, mix_mus = mix_mus, mix_sigmas = mix_sigmas,
    mix_skews = mix_skews, mix_skurts = mix_skurts, mix_fifths = mix_fifths,
    mix_sixths = mix_sixths, mix_Six = mix_Six, marginal = marginal,
    support = support, lam = lam, p_zip = p_zip[2],
    rho = Rey[-c(8:9), -c(8:9)], seed = seed,
    epsilon = 0.01, errorloop = TRUE)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})

means <- c(L[1], L[1])
vars <- c(L[2]^2, L[2]^2)
skews <- c(L[3], L[3])
skurts <- c(L[4], L[4])
fifths <- c(L[5], L[5])
sixths <- c(L[6], L[6])
Six <- list(1.75, 1.75)

test_that("works for Fleishman method: 1 ordinal, 2 iid continuous,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 2, k_mix = 0, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts,
    marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey[-c(4:5), -c(4:5)], seed = seed,
    epsilon = 0.01)$constants[1, "c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 2, k_mix = 0, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts,
    marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey[-c(4:5), -c(4:5)], seed = seed, epsilon = 0.01,
    cstart = cstartF)$constants[1, "c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 2, k_mix = 0, k_pois, k_nb,
    "Fleishman", means, vars, skews, skurts,
    marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey[-c(4:5), -c(4:5)], seed = seed, epsilon = 0.01,
    errorloop = TRUE)$constants[1, "c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 2 iid continuous,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 2, k_mix = 0, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey[-c(4:5), -c(4:5)], seed = seed,
    epsilon = 0.01)$constants[1, "c5"], 0.0000006124845, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 2, k_mix = 0, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey[-c(4:5), -c(4:5)], seed = seed, epsilon = 0.01,
    cstart = cstartP)$constants[1, "c5"], 0.0000006124845, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 2, k_mix = 0, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey[-c(4:5), -c(4:5)], seed = seed, epsilon = 0.01,
    errorloop = TRUE)$constants[1, "c5"], 0.0000006124845, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

mix_pis <- list(c(0.3, 0.2, 0.5))
mix_mus <- list(c(L[1], L[1], B[1]))
mix_sigmas <- list(c(L[2], L[2], B[2]))
mix_skews <- list(c(L[3], L[3], B[3]))
mix_skurts <- list(c(L[4], L[4], B[4]))
mix_fifths <- list(c(L[5], L[5], B[5]))
mix_sixths <- list(c(L[6], L[6], B[6]))
mix_Six <- list(1.75, 1.75, 0.03)
Mstcum <- calc_mixmoments(mix_pis[[1]], mix_mus[[1]], mix_sigmas[[1]],
  mix_skews[[1]], mix_skurts[[1]], mix_fifths[[1]], mix_sixths[[1]])

means <- Mstcum[1]
vars <- Mstcum[2]^2

test_that("works for Fleishman method: 1 ordinal, 1 with 2 iid mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey[-2, -2], seed = seed, epsilon = 0.01)$constants[1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey[-2, -2], seed = seed, epsilon = 0.01,
    cstart = cstartF)$constants[1, "c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Fleishman", means, vars, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey[-2, -2], seed = seed, epsilon = 0.01,
    errorloop = TRUE)$constants[1, "c3"], 0.03605955, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method: 1 ordinal, 1 with 2 iid mixture,
          1 Poisson, 1 ZIP, 1 NB, 1 ZINB", {
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, mix_fifths = mix_fifths, mix_sixths = mix_sixths,
    mix_Six = mix_Six, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, prob = prob, p_zinb = p_zinb[2],
    rho = Rey[-2, -2], seed = seed, epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, mix_fifths = mix_fifths, mix_sixths = mix_sixths,
    mix_Six = mix_Six, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip, size = size, mu = mu, p_zinb = p_zinb,
    rho = Rey[-2, -2], seed = seed, epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(corrvar2(n, k_cat, k_cont = 0, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, mix_pis = mix_pis,
    mix_mus = mix_mus, mix_sigmas = mix_sigmas, mix_skews = mix_skews,
    mix_skurts = mix_skurts, mix_fifths = mix_fifths, mix_sixths = mix_sixths,
    mix_Six = mix_Six, marginal = marginal, support = support,
    lam = lam, p_zip = p_zip[2], size = size, mu = mu, p_zinb = p_zinb[2],
    rho = Rey[-2, -2], seed = seed, epsilon = 0.01)$constants[1, "c5"],
    0.0000006124845, tolerance = tol, check.attributes = FALSE), TRUE)
})
