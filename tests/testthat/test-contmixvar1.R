skip_on_cran()
library("SimCorrMix")
context("Simulate one continuous mixture variable")

options(scipen = 999)
tol <- 1e-5
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
mix_Six <- list(seq(0.01, 10, 0.01), c(0.01, 0.02, 0.03),
  seq(0.01, 10, 0.01))
Bstcum <- calc_mixmoments(mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths)

n <- 25
set.seed(1234)
cstart1 <- runif(n, min = -2, max = 2)
cstart2 <- runif(n, min = -1, max = 1)
cstart3 <- runif(n, min = -0.5, max = 0.5)
cstartF <- cbind(cstart1, cstart2, cstart3)
cstartF <- list(cstartF, cstartF, cstartF)

set.seed(1234)
cstart1 <- runif(n, min = -2, max = 2)
cstart2 <- runif(n, min = -1, max = 1)
cstart3 <- runif(n, min = -1, max = 1)
cstart4 <- runif(n, min = -0.025, max = 0.025)
cstart5 <- runif(n, min = -0.025, max = 0.025)
cstartP <- cbind(cstart1, cstart2, cstart3, cstart4, cstart5)
cstartP <- list(cstartP, cstartP, cstartP)

test_that("works for Fleishman method", {
  expect_equal(all.equal(contmixvar1(n = 10000, "Fleishman", Bstcum[1],
    Bstcum[2]^2, mix_pis, mix_mus, mix_sigmas, mix_skews,
    mix_skurts)$constants[1, "c3"], -0.026275030, tolerance = tol,
    check.attributes = FALSE), TRUE)
})

test_that("works for Fleishman method given cstart", {
  expect_equal(all.equal(contmixvar1(n = 10000, "Fleishman", Bstcum[1],
    Bstcum[2]^2, mix_pis, mix_mus, mix_sigmas, mix_skews,
    mix_skurts, cstart = cstartF)$constants[1, "c3"], -0.026275030,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method", {
  expect_equal(all.equal(contmixvar1(n = 10000, "Polynomial", Bstcum[1],
    Bstcum[2]^2, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
    mix_fifths, mix_sixths)$constants[1, "c5"], 0.00061977431,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method with a Six list", {
  expect_equal(all.equal(contmixvar1(n = 10000, "Polynomial", Bstcum[1],
    Bstcum[2]^2, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
    mix_fifths, mix_sixths, mix_Six)$constants[1, "c5"], 0.00061977431,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("works for Polynomial method given cstart", {
  expect_equal(all.equal(contmixvar1(n = 10000, "Polynomial", Bstcum[1],
    Bstcum[2]^2, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
    mix_fifths, mix_sixths, mix_Six, cstart = cstartP)$constants[1, "c5"],
    0.00061977431, tolerance = tol, check.attributes = FALSE), TRUE)
})
