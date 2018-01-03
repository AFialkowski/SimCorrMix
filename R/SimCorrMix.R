#' @title Simulation of Correlated Data of Multiple Variable Types including Continuous and Count Mixture Distributions
#'
#' @description \pkg{SimCorrMix} generates continuous (normal, non-normal, or mixture distributions), binary, ordinal, and count
#'     (Poisson or Negative Binomial, regular or zero-inflated) variables with a specified correlation matrix, or one continuous variable
#'     with a mixture distribution.  This package can be used to simulate data sets that mimic real-world clinical or genetic data sets
#'     (i.e. plasmodes, as in Vaughan et al., 2009, \doi{10.1016/j.csda.2008.02.032}).  The methods extend those found in the
#'     \pkg{SimMultiCorrData} package.  Standard normal variables with an imposed intermediate correlation matrix are transformed to
#'     generate the desired distributions.  Continuous variables are simulated using either Fleishman's third-order
#'     (\doi{10.1007/BF02293811}) or Headrick's fifth-order (\doi{10.1016/S0167-9473(02)00072-5}) power method transformation (PMT).
#'     Non-mixture distributions require the user to specify mean, variance, skewness, standardized kurtosis, and standardized fifth and
#'     sixth cumulants.  Mixture distributions require these inputs for the component distributions plus the mixing probabilities.  Simulation
#'     occurs at the component-level for continuous mixture distributions.  The target correlation matrix is specified in terms of
#'     correlations with components of continuous mixture variables.  However, the package provides functions to approximate expected
#'     correlations with continuous mixture variables given target correlations with the components.  Binary and ordinal variables are simulated using a modification of
#'     \code{\link[GenOrd]{GenOrd-package}}'s \code{\link[GenOrd]{ordsample}} function.  Count variables are simulated using the inverse
#'     CDF method.  There are two simulation pathways which calculate intermediate correlations involving count variables differently.
#'     Correlation Method 1 adapts Yahav and Shmueli's 2012 method (\doi{10.1002/asmb.901}).  Correlation Method 2 adapts Barbiero and
#'     Ferrari's 2015 modification of \code{\link[GenOrd]{GenOrd-package}} (\doi{10.1002/asmb.2072}).  The optional error loop may be used
#'     to improve the accuracy of the final correlation matrix.  The package also provides functions to calculate the standardized
#'     cumulants of continuous mixture distributions, check parameter inputs, calculate feasible correlation boundaries, and summarize and plot simulated variables.
#'
#'
#' @seealso Useful link: \url{https://github.com/AFialkowski/SimCorrMix}
#' @section Vignettes:
#' There are several vignettes which accompany this package to help the user understand the simulation and analysis methods.
#'
#' 1) \bold{Calculation of Correlation Boundaries} explains how the feasible correlation boundaries are calculated for each of the two
#' simulation pathways.
#'
#' 2) \bold{Comparison of Correlation Methods 1 and 2} describes the two simulation pathways that can be followed for generation of
#' correlated data.
#'
#' 3) \bold{Continuous Mixture Distributions} demonstrates how to simulate one continuous mixture variable using
#' \code{\link[SimCorrMix]{contmixvar1}} and gives a step-by-step guideline for comparing a simulated distribution to the target
#' distribution.
#'
#' 4) \bold{Error Loop Algorithm} details the algorithm involved in the optional error loop that helps to minimize correlation errors.
#'
#' 5) \bold{Expected Cumulants and Correlations for Continuous Mixture Variables} derives the equations used by the function
#' \code{\link[SimCorrMix]{calc_mixmoments}} to find the mean, standard deviation, skew, standardized kurtosis, and standardized fifth
#' and sixth cumulants for a continuous mixture variable.  The vignette also explains how the functions
#' \code{\link[SimCorrMix]{rho_M1M2}} and \code{\link[SimCorrMix]{rho_M1Y}} approximate the expected correlations with continuous mixture
#' variables based on the target correlations with the components.
#'
#' 6) \bold{Overall Workflow for Generation of Correlated Data} gives a step-by-step guideline to follow with an example containing
#' continuous non-mixture and mixture, ordinal, zero-inflated Poisson, and zero-inflated Negative Binomial variables.  It executes both
#' correlated data simulation functions with and without the error loop.
#'
#' 7) \bold{Variable Types} describes the different types of variables that can be simulated in \pkg{SimCorrMix}.
#'
#'
#' @section Functions:
#' This package contains 3 \emph{simulation} functions:
#'
#' \code{\link[SimCorrMix]{contmixvar1}}, \code{\link[SimCorrMix]{corrvar}}, and \code{\link[SimCorrMix]{corrvar2}}
#'
#' 4 data description (\emph{summary}) function:
#'
#' \code{\link[SimCorrMix]{calc_mixmoments}}, \code{\link[SimCorrMix]{summary_var}}, \code{\link[SimCorrMix]{rho_M1M2}}, \code{\link[SimCorrMix]{rho_M1Y}}
#'
#' 2 \emph{graphing} functions:
#'
#' \code{\link[SimCorrMix]{plot_simpdf_theory}}, \code{\link[SimCorrMix]{plot_simtheory}}
#'
#' 3 \emph{support} functions:
#'
#' \code{\link[SimCorrMix]{validpar}}, \code{\link[SimCorrMix]{validcorr}}, \code{\link[SimCorrMix]{validcorr2}}
#'
#' and 16 \emph{auxiliary} functions (should not normally be called by the user, but are called by other functions):
#'
#' \code{\link[SimCorrMix]{corr_error}}, \code{\link[SimCorrMix]{intercorr}}, \code{\link[SimCorrMix]{intercorr2}},
#' \code{\link[SimCorrMix]{intercorr_cat_nb}}, \code{\link[SimCorrMix]{intercorr_cat_pois}}, \cr
#' \code{\link[SimCorrMix]{intercorr_cont_nb}}, \code{\link[SimCorrMix]{intercorr_cont_nb2}},
#' \code{\link[SimCorrMix]{intercorr_cont_pois}}, \code{\link[SimCorrMix]{intercorr_cont_pois2}}, \cr
#' \code{\link[SimCorrMix]{intercorr_cont}}, \code{\link[SimCorrMix]{intercorr_nb}}, \code{\link[SimCorrMix]{intercorr_pois}},
#' \code{\link[SimCorrMix]{intercorr_pois_nb}}, \code{\link[SimCorrMix]{maxcount_support}},
#' \code{\link[SimCorrMix]{ord_norm}}, \code{\link[SimCorrMix]{norm_ord}}
#'
#' @docType package
#' @name SimCorrMix
#' @references
#' Amatya A & Demirtas H (2015). Simultaneous generation of multivariate mixed data with Poisson and normal marginals.
#'     Journal of Statistical Computation and Simulation, 85(15):3129-39. \doi{10.1080/00949655.2014.953534}.
#'
#' Barbiero A & Ferrari PA (2015). Simulation of correlated Poisson variables. Applied Stochastic Models in
#'     Business and Industry, 31:669-80. \doi{10.1002/asmb.2072}.
#'
#' Barbiero A & Ferrari PA (2015). GenOrd: Simulation of Discrete Random Variables with Given
#'     Correlation Matrix and Marginal Distributions. R package version 1.4.0. \cr \url{https://CRAN.R-project.org/package=GenOrd}
#'
#' Berend H (2017). nleqslv: Solve Systems of Nonlinear Equations. R package version 3.2.
#'     \url{https://CRAN.R-project.org/package=nleqslv}
#'
#' Carnell R (2017). triangle: Provides the Standard Distribution Functions for the Triangle Distribution. R package version 0.11.
#'     \url{https://CRAN.R-project.org/package=triangle}.
#'
#' Davenport JW, Bezder JC, & Hathaway RJ (1988). Parameter Estimation for Finite Mixture Distributions.
#'     Computers & Mathematics with Applications, 15(10):819-28.
#'
#' Demirtas H (2006). A method for multivariate ordinal data generation given marginal distributions and correlations. Journal of Statistical
#'     Computation and Simulation, 76(11):1017-1025. \cr \doi{10.1080/10629360600569246}.
#'
#' Demirtas H (2014). Joint Generation of Binary and Nonnormal Continuous Data. Biometrics & Biostatistics, S12.
#'
#' Demirtas H & Hedeker D (2011). A practical way for computing approximate lower and upper correlation bounds.
#'     American Statistician, 65(2):104-109. \doi{10.1198/tast.2011.10090}.
#'
#' Demirtas H, Hedeker D, & Mermelstein RJ (2012). Simulation of massive public health data by power polynomials.
#'     Statistics in Medicine, 31(27):3337-3346. \doi{10.1002/sim.5362}.
#'
#' Emrich LJ & Piedmonte MR (1991). A Method for Generating High-Dimensional Multivariate Binary Variables. The American Statistician, 45(4): 302-4.
#'     \doi{10.1080/00031305.1991.10475828}.
#'
#' Everitt BS (1996). An Introduction to Finite Mixture Distributions. Statistical Methods in Medical Research, 5(2):107-127. \doi{10.1177/096228029600500202}.
#'
#' Ferrari PA & Barbiero A (2012). Simulating ordinal data. Multivariate Behavioral Research, 47(4): 566-589.
#'     \doi{10.1080/00273171.2012.692630}.
#'
#' Fialkowski AC (2017). SimMultiCorrData: Simulation of Correlated Data with Multiple Variable Types. R package version 0.2.1.
#'     \url{https://CRAN.R-project.org/package=SimMultiCorrData}.
#'
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43:521-532. \doi{10.1007/BF02293811}.
#'
#' Frechet M (1951). Sur les tableaux de correlation dont les marges sont donnees.  Ann. l'Univ. Lyon SectA, 14:53-77.
#'
#' Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#'     Non-normal Distributions. Computational Statistics & Data Analysis, 40(4):685-711. \doi{10.1016/S0167-9473(02)00072-5}.
#'     (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3(1):65-71. \doi{10.22237/jmasm/1083370080}.
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77:229-249. \doi{10.1080/10629360600605065}.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64:25-35. \doi{10.1007/BF02294317}.
#'
#' Headrick TC, Sawilowsky SS (2002). Weighted Simplex Procedures for Determining Boundary Points and Constants for the
#'     Univariate and Multivariate Power Methods. Journal of Educational and Behavioral Statistics, 25:417-436. \doi{10.3102/10769986025004417}.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3):1 - 17. \cr \doi{10.18637/jss.v019.i03}.
#'
#' Higham N (2002). Computing the nearest correlation matrix - a problem from finance; IMA Journal of Numerical Analysis 22:329-343.
#'
#' Hoeffding W. Scale-invariant correlation theory. In: Fisher NI, Sen PK, editors. The collected works of Wassily Hoeffding.
#'     New York: Springer-Verlag; 1994. p. 57-107.
#'
#' Ismail N & Zamani H (2013). Estimation of Claim Count Data Using Negative Binomial, Generalized Poisson, Zero-Inflated Negative Binomial and
#'     Zero-Inflated Generalized Poisson Regression Models. Casualty Actuarial Society E-Forum 41(20):1-28.
#'
#' Lambert D (1992). Zero-Inflated Poisson Regression, with an Application to Defects in Manufacturing. Technometrics 34(1):1-14.
#'
#' Olsson U, Drasgow F, & Dorans NJ (1982). The Polyserial Correlation Coefficient. Psychometrika, 47(3):337-47.
#'     \doi{10.1007/BF02294164}.
#'
#' Pearson RK (2011). Exploring Data in Engineering, the Sciences, and Medicine. In. New York: Oxford University Press.
#'
#' Schork NJ, Allison DB, & Thiel B (1996). Mixture Distributions in Human Genetics Research. Statistical Methods in Medical Research,
#'     5:155-178. \doi{10.1177/096228029600500204}.
#'
#' Vale CD & Maurelli VA (1983). Simulating Multivariate Nonnormal Distributions. Psychometrika, 48:465-471. \doi{10.1007/BF02293687}.
#'
#' Varadhan R, Gilbert PD (2009). BB: An R Package for Solving a Large System of Nonlinear Equations and for
#'     Optimizing a High-Dimensional Nonlinear Objective Function, J. Statistical Software, 32(4). \doi{10.18637/jss.v032.i04}.
#'     \url{http://www.jstatsoft.org/v32/i04/}
#'
#' Vaughan LK, Divers J, Padilla M, Redden DT, Tiwari HK, Pomp D, Allison DB (2009). The use of plasmodes as a supplement to simulations:
#'     A simple example evaluating individual admixture estimation methodologies. Comput Stat Data Anal, 53(5):1755-66.
#'     \doi{10.1016/j.csda.2008.02.032}.
#'
#' Yahav I & Shmueli G (2012). On Generating Multivariate Poisson Data in Management Science Applications. Applied Stochastic
#'     Models in Business and Industry, 28(1):91-102. \doi{10.1002/asmb.901}.
#'
#' Yee TW (2017). VGAM: Vector Generalized Linear and Additive Models. \cr \url{https://CRAN.R-project.org/package=VGAM}.
#'
#' Zhang X, Mallick H, & Yi N (2016). Zero-Inflated Negative Binomial Regression for Differential Abundance Testing in Microbiome
#'     Studies. Journal of Bioinformatics and Genomics 2(2):1-9. \doi{10.18454/jbg.2016.2.2.1}.
#'
NULL
