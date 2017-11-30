## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, message=FALSE, fig.width = 6, fig.height = 4, cache = FALSE)

## ---- include=FALSE------------------------------------------------------
library("SimCorrMix")
library("printr")

## ------------------------------------------------------------------------
options(scipen = 999)
seed <- 276
n <- 10000

# Continuous variables
L <- calc_theory("Logistic", c(0, 1))
C <- calc_theory("Chisq", 4)
B <- calc_theory("Beta", c(4, 1.5))

# Non-mixture variables
skews <- rep(L[3], 2)
skurts <- rep(L[4], 2)
fifths <- rep(L[5], 2)
sixths <- rep(L[6], 2)
Six <- list(1.75, 1.75)

# Mixture variables
mix_pis <- list(c(0.4, 0.6), c(0.3, 0.2, 0.5))
mix_mus <- list(c(-2, 2), c(L[1], C[1], B[1]))
mix_sigmas <- list(c(1, 1), c(L[2], C[2], B[2]))
mix_skews <- list(rep(0, 2), c(L[3], C[3], B[3]))
mix_skurts <- list(rep(0, 2), c(L[4], C[4], B[4]))
mix_fifths <- list(rep(0, 2), c(L[5], C[5], B[5]))
mix_sixths <- list(rep(0, 2), c(L[6], C[6], B[6]))
mix_Six <- list(list(NULL, NULL), list(1.75, NULL, 0.03))
Nstcum <- calc_mixmoments(mix_pis[[1]], mix_mus[[1]], mix_sigmas[[1]], 
  mix_skews[[1]], mix_skurts[[1]], mix_fifths[[1]], mix_sixths[[1]])
Mstcum <- calc_mixmoments(mix_pis[[2]], mix_mus[[2]], mix_sigmas[[2]], 
  mix_skews[[2]], mix_skurts[[2]], mix_fifths[[2]], mix_sixths[[2]])

means <- c(L[1], L[1], Nstcum[1], Mstcum[1])
vars <- c(L[2]^2, L[2]^2, Nstcum[2]^2, Mstcum[2]^2)

marginal <- list(0.3)
support <- list(c(0, 1))
lam <- 0.5
p_zip <- 0.1
size <- 2
prob <- 0.75
mu <- size * (1 - prob)/prob
p_zinb <- 0.2

k_cat <- length(marginal) 
k_cont <- length(Six)
k_mix <- length(mix_pis)
k_comp <- sum(unlist(lapply(mix_pis, length)))
k_pois <- length(lam)
k_nb <- length(size)
k_total <- k_cat + k_cont + k_comp + k_pois + k_nb

Rey <- matrix(0.35, k_total, k_total)
diag(Rey) <- 1
rownames(Rey) <- colnames(Rey) <- c("O1", "C1", "C2", "M1_1", "M1_2", "M2_1", 
  "M2_2", "M2_3", "P1", "NB1")
Rey["M1_1", "M1_2"] <- Rey["M1_2", "M1_1"] <- 0
Rey["M2_1", "M2_2"] <- Rey["M2_2", "M2_1"] <- Rey["M2_1", "M2_3"] <- 
  Rey["M2_3", "M2_1"] <- Rey["M2_2", "M2_3"] <- Rey["M2_3", "M2_2"] <- 0

## ------------------------------------------------------------------------
validpar(k_cat, k_cont, k_mix, k_pois, k_nb, "Polynomial", 
  means, vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, 
  mix_sigmas, mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, 
  marginal, lam, p_zip, size, prob, mu = NULL, p_zinb, rho = Rey)

## ------------------------------------------------------------------------
Lower_third <- calc_lower_skurt(method = "Fleishman", skews = C[3], 
  Skurt = seq(1.161, 1.17, 0.001), seed = 104)
knitr::kable(Lower_third$Min[, c("skew", "valid.pdf", "skurtosis")], 
  row.names = FALSE, caption = "Third-Order Lower Skurtosis Bound")

## ------------------------------------------------------------------------
Lower_fifth <- calc_lower_skurt(method = "Polynomial", skews = C[3], 
  fifths = C[5], sixths = C[6], Skurt = seq(0.022, 0.03, 0.001), seed = 104)
knitr::kable(Lower_fifth$Min[, c("skew", "fifth", "sixth", "valid.pdf", 
  "skurtosis")], row.names = FALSE, 
  caption = "Fifth-Order Lower Skurtosis Bound")

## ------------------------------------------------------------------------
valid1 <- validcorr(n, k_cat, k_cont, k_mix, k_pois, k_nb, "Polynomial", 
  means, vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, 
  mix_sigmas, mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, 
  marginal, lam, p_zip, size, prob, mu = NULL, p_zinb, Rey, seed)
valid1$valid.rho

## ------------------------------------------------------------------------
Sim1 <- corrvar(n, k_cat, k_cont, k_mix, k_pois, k_nb, 
  "Polynomial", means, vars, skews, skurts, fifths, sixths, Six, 
  mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, mix_fifths, 
  mix_sixths, mix_Six, marginal, support, lam, p_zip, size, prob, 
  mu = NULL, p_zinb, Rey, seed, epsilon = 0.01)
Sim1_error <- abs(Rey - Sim1$rho_calc)

## ------------------------------------------------------------------------
summary(as.numeric(Sim1_error))

## ------------------------------------------------------------------------
rho_mix <- Sim1$rho_mix
rownames(rho_mix) <- c("01", "C1", "C2", "M1", "M2", "P1", "NB1")
colnames(rho_mix) <- rownames(rho_mix)
rho_mix

## ------------------------------------------------------------------------
p_M11M21 <- p_M11M22 <- p_M11M23 <- 0.35
p_M12M21 <- p_M12M22 <- p_M12M23 <- 0.35
p_M1M2 <- matrix(c(p_M11M21, p_M11M22, p_M11M23, p_M12M21, p_M12M22, p_M12M23), 
  2, 3, byrow = TRUE)
rhoM1M2 <- rho_M1M2(mix_pis, mix_mus, mix_sigmas, p_M1M2)

## ------------------------------------------------------------------------
p_M11C1 <- p_M12C1 <- 0.35
p_M1C1 <- c(p_M11C1, p_M12C1)
rho_M1C1 <- rho_M1Y(mix_pis[[1]], mix_mus[[1]], mix_sigmas[[1]], p_M1C1)

## ------------------------------------------------------------------------
p_M21C1 <- p_M22C1 <- p_M23C1 <- 0.35
p_M2C1 <- c(p_M21C1, p_M22C1, p_M23C1)
rho_M2C1 <- rho_M1Y(mix_pis[[2]], mix_mus[[2]], mix_sigmas[[2]], p_M2C1)

## ------------------------------------------------------------------------
Sim1_EL <- corrvar(n, k_cat, k_cont, k_mix, k_pois, k_nb, 
  "Polynomial", means, vars, skews, skurts, fifths, sixths, Six, 
  mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, mix_fifths, 
  mix_sixths, mix_Six, marginal, support, lam, p_zip, size, prob, 
  mu = NULL, p_zinb, Rey, seed, errorloop = TRUE, epsilon = 0.01)
EL1_error <- abs(Rey - Sim1_EL$rho_calc)

## ------------------------------------------------------------------------
summary(as.numeric(EL1_error))

## ------------------------------------------------------------------------
rho_mixEL <- Sim1_EL$rho_mix
rownames(rho_mixEL) <- c("01", "C1", "C2", "M1", "M2", "P1", "NB1")
colnames(rho_mixEL) <- rownames(rho_mixEL)
rho_mixEL

## ------------------------------------------------------------------------
Sim1_EL$valid.pdf
Sim1_EL$sixth_correction

## ------------------------------------------------------------------------
target_sum <- Sim1_EL$target_sum
cont_sum <- Sim1_EL$cont_sum
rownames(target_sum) <- rownames(cont_sum) <- c("C1", "C2", "M1_1", "M1_2", 
  "M2_1", "M2_2", "M2_3")
knitr::kable(target_sum, digits = 5, row.names = TRUE, 
  caption = "Summary of Target Distributions")
knitr::kable(cont_sum[, -c(2, 5:7)], digits = 5, row.names = TRUE, 
  caption = "Summary of Simulated Distributions")

## ------------------------------------------------------------------------
target_mix <- Sim1_EL$target_mix
mix_sum <- Sim1_EL$mix_sum
rownames(target_mix) <- rownames(mix_sum) <- c("M1", "M2")
knitr::kable(target_mix, digits = 5, row.names = TRUE, 
  caption = "Summary of Target Distributions")
knitr::kable(mix_sum[, -c(2, 5:7)], digits = 5, row.names = TRUE, 
  caption = "Summary of Simulated Distributions")

## ------------------------------------------------------------------------
Nplot <- plot_simpdf_theory(sim_y = Sim1_EL$Y_mix[, 1], ylower = -10, 
  yupper = 10, title = "PDF of Mixture of N(-2, 1) and N(2, 1) Distributions",
  fx = function(x) mix_pis[[1]][1] * dnorm(x, mix_mus[[1]][1], 
    mix_sigmas[[1]][1]) + mix_pis[[1]][2] * dnorm(x, mix_mus[[1]][2], 
    mix_sigmas[[1]][2]), lower = -Inf, upper = Inf)
Nplot
Mplot <- plot_simpdf_theory(sim_y = Sim1_EL$Y_mix[, 2], 
  title = paste("PDF of Mixture of Logistic(0, 1), Chisq(4),", 
    "\nand Beta(4, 1.5) Distributions", sep = ""),
  fx = function(x) mix_pis[[2]][1] * dlogis(x, 0, 1) + mix_pis[[2]][2] * 
    dchisq(x, 4) + mix_pis[[2]][3] * dbeta(x, 4, 1.5), 
  lower = -Inf, upper = Inf)
Mplot

## ------------------------------------------------------------------------
knitr::kable(Sim1_EL$ord_sum, caption = "Summary of Ordinal Variables")
knitr::kable(Sim1_EL$pois_sum[, -c(2, 9:11)], 
  caption = "Summary of Poisson Variables")
Pplot <- plot_simpdf_theory(sim_y = Sim1_EL$Y_pois[, 1], 
  title = "PMF of Zero-Inflated Poisson Distribution", Dist = "Poisson", 
  params = c(lam, p_zip), cont_var = FALSE)
Pplot

## ------------------------------------------------------------------------
knitr::kable(Sim1_EL$nb_sum[, -c(2, 10:12)], 
  caption = "Summary of Negative Binomial Variables")
NBplot <- plot_simtheory(sim_y = Sim1_EL$Y_nb[, 1], 
  title = "Simulated Zero-Inflated NB Values", 
  Dist = "Negative_Binomial", params = c(size, mu, p_zinb), 
  cont_var = FALSE)
NBplot

## ------------------------------------------------------------------------
pois_eps <- 0.0001
nb_eps <- 0.0001
valid2 <- validcorr2(n, k_cat, k_cont, k_mix, k_pois, k_nb, "Polynomial", 
  means, vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, 
  mix_sigmas, mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, 
  lam, p_zip, size, prob = NULL, mu, p_zinb, pois_eps, nb_eps, Rey, seed)
valid2$valid.rho

## ------------------------------------------------------------------------
Sim2 <- corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb, 
  "Polynomial", means, vars, skews, skurts, fifths, sixths, Six, 
  mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, mix_fifths, 
  mix_sixths, mix_Six, marginal, support, lam, p_zip, size, prob = NULL, mu, 
  p_zinb, pois_eps, nb_eps, Rey, seed, epsilon = 0.01)
Sim2_error <- abs(Rey - Sim2$rho_calc)

## ------------------------------------------------------------------------
summary(as.numeric(Sim2_error))

## ------------------------------------------------------------------------
rho_mix <- Sim2$rho_mix
rownames(rho_mix) <- c("01", "C1", "C2", "M1", "M2", "P1", "NB1")
colnames(rho_mix) <- rownames(rho_mix)
rho_mix

## ------------------------------------------------------------------------
Sim2_EL <- corrvar2(n, k_cat, k_cont, k_mix, k_pois, k_nb, 
  "Polynomial", means, vars, skews, skurts, fifths, sixths, Six, mix_pis, 
  mix_mus, mix_sigmas, mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, 
  marginal, support, lam, p_zip, size, prob = NULL, mu, p_zinb, 
  pois_eps, nb_eps, Rey, seed, errorloop = TRUE, epsilon = 0.01)
EL2_error <- abs(Rey - Sim2_EL$rho_calc)

## ------------------------------------------------------------------------
summary(as.numeric(EL2_error))

## ------------------------------------------------------------------------
rho_mixEL <- Sim2_EL$rho_mix
rownames(rho_mixEL) <- c("01", "C1", "C2", "M1", "M2", "P1", "NB1")
colnames(rho_mixEL) <- rownames(rho_mixEL)
rho_mixEL

## ------------------------------------------------------------------------
Sim1_EL$valid.pdf
Sim1_EL$sixth_correction

## ------------------------------------------------------------------------
target_sum <- Sim2$target_sum
cont_sum <- Sim2$cont_sum
rownames(target_sum) <- rownames(cont_sum) <- c("C1", "C2", "M1_1", "M1_2", 
  "M2_1", "M2_2", "M2_3")
knitr::kable(target_sum, digits = 5, row.names = TRUE, 
  caption = "Summary of Target Distributions")
knitr::kable(cont_sum[, -c(2, 5:7)], digits = 5, row.names = TRUE, 
  caption = "Summary of Simulated Distributions")

## ------------------------------------------------------------------------
target_mix <- Sim2$target_mix
mix_sum <- Sim2$mix_sum
rownames(target_mix) <- rownames(mix_sum) <- c("M1", "M2")

knitr::kable(target_mix, digits = 5, row.names = TRUE, 
             caption = "Summary of Target Distributions")
knitr::kable(mix_sum[, -c(2, 5:7)], digits = 5, row.names = TRUE, 
             caption = "Summary of Simulated Distributions")

## ------------------------------------------------------------------------
Nplot <- plot_simpdf_theory(sim_y = Sim2$Y_mix[, 1], ylower = -10, 
  yupper = 10, title = "Mixture of N(-2, 1) and N(2, 1) Distributions",
  fx = function(x) mix_pis[[1]][1] * dnorm(x, mix_mus[[1]][1], 
    mix_sigmas[[1]][1]) + mix_pis[[1]][2] * dnorm(x, mix_mus[[1]][2], 
    mix_sigmas[[1]][2]), lower = -Inf, upper = Inf)
Nplot
Mplot <- plot_simpdf_theory(sim_y = Sim2$Y_mix[, 2], 
  title = paste("Mixture of Logistic(0, 1), Chisq(4),", 
    "\nand Beta(4, 1.5) Distributions", sep = ""),
  fx = function(x) mix_pis[[2]][1] * dlogis(x, 0, 1) + mix_pis[[2]][2] * 
    dchisq(x, 4) + mix_pis[[2]][3] * dbeta(x, 4, 1.5), 
  lower = -Inf, upper = Inf)
Mplot

## ------------------------------------------------------------------------
knitr::kable(Sim2$ord_sum, caption = "Summary of Ordinal Variables")
knitr::kable(Sim2$pois_sum[, -c(2, 9:11)], 
  caption = "Summary of Poisson Variables")
Pplot <- plot_simpdf_theory(sim_y = Sim2$Y_pois[, 1], 
  title = "PMF of Zero-Inflated Poisson Distribution", Dist = "Poisson", 
  params = c(lam, p_zip), cont_var = FALSE)
Pplot

## ------------------------------------------------------------------------
knitr::kable(Sim2$nb_sum[, -c(2, 10:12)], 
  caption = "Summary of Negative Binomial Variables")
NBplot <- plot_simtheory(sim_y = Sim2$Y_nb[, 1], 
  title = "Simulated Zero-Inflated NB Values", 
  Dist = "Negative_Binomial", params = c(size, mu, p_zinb), 
  cont_var = FALSE)
NBplot

