skip_on_cran()
library("SimCorrMix")
context("Test summary functions")

tol <- 1e-5
L <- calc_theory("Logistic", c(0, 1))
C <- calc_theory("Chisq", 4)
B <- calc_theory("Beta", c(4, 1.5))
mix_pis <- list(c(0.4, 0.6), c(0.3, 0.2, 0.5))
mix_mus <- list(c(-2, 2), c(L[1], C[1], B[1]))
mix_sigmas <- list(c(1, 1), c(L[2], C[2], B[2]))
mix_skews <- list(rep(0, 2), c(L[3], C[3], B[3]))
mix_skurts <- list(rep(0, 2), c(L[4], C[4], B[4]))
mix_fifths <- list(rep(0, 2), c(L[5], C[5], B[5]))
mix_sixths <- list(rep(0, 2), c(L[6], C[6], B[6]))
p_M11M21 <- p_M11M22 <- p_M11M23 <- 0.35
p_M12M21 <- p_M12M22 <- p_M12M23 <- 0.35
p_M1M2 <- matrix(c(p_M11M21, p_M11M22, p_M11M23, p_M12M21, p_M12M22, p_M12M23),
  2, 3, byrow = TRUE)
p_M11C1 <- p_M12C1 <- 0.35
p_M1C1 <- c(p_M11C1, p_M12C1)

test_that("calculate mixture cumulants for Fleishman", {
  expect_equal(all.equal(calc_mixmoments(mix_pis[[1]], mix_mus[[1]],
    mix_sigmas[[1]], mix_skews[[1]], mix_skurts[[1]])[2], 2.20,
    tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("calculate mixture cumulants for Polynomial", {
  expect_equal(all.equal(calc_mixmoments(mix_pis[[1]], mix_mus[[1]],
    mix_sigmas[[1]], mix_skews[[1]], mix_skurts[[1]], mix_fifths[[1]],
    mix_sixths[[1]])[2], 2.20, tolerance = tol, check.attributes = FALSE),
    TRUE)
})

test_that("calculate correlation between M1 and M2", {
  expect_equal(all.equal(rho_M1M2(mix_pis, mix_mus, mix_sigmas, p_M1M2),
    0.08773416, tolerance = tol, check.attributes = FALSE), TRUE)
})

test_that("calculate correlation between M1 and Y", {
  expect_equal(all.equal(rho_M1Y(mix_pis[[1]], mix_mus[[1]], mix_sigmas[[1]],
    p_M1C1), 0.1590909, tolerance = tol, check.attributes = FALSE), TRUE)
})
