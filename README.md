# Nonparametric-L-moment-homogeneity-test
***
R code implementing the **Nonparametric L-moment homogeneity test** proposed in the article:
  
> Masselot, P., Chebana, F., Ouarda, T.B.M.J., 2017. Fast and direct nonparametric procedures in the L-moment homogeneity test. Stochastic Environmental Research and Risk Assessment 31, 509–522. [https://doi.org/10.1007/s00477-016-1248-0](https://doi.org/10.1007/s00477-016-1248-0)

This methodological paper proposes a nonparametric extension of the Hosking-Wallis homogeneity test, extensively used in regional frequency analysis. This procedure is used to test whether a set of hydrological sites (a region) can be considered as having the same distribution up to a scale factor.

In the Hosking-Wallis test, the null distribution is simulated from a flexible four parameters Kappa distribution. This results in an important amount of uncertainty due to the four parameter estimation, performed through numerical optimization for which convergence is not guaranteed. 

In the proposed version, the test statistic remain unchanged, but its null distribution is now simulated using a nonparametric procedure. This has the advantage of avoiding the estimation of a parametric distribution. Four different nonparametric procedures are available in the function to simulate homogeneous regions:
* *Permutations*: permutes the observations from the sites;
* *Bootstrap*: pools the observations of all sites and then create new sites by drawing observations with replacement;
* *Polya*: Similar to bootstrap, but each time an observation is drawn its value is duplicated in the pooled sample. Two variants of this scheme are available:
  + *Site-wise*: each simulated site is drawn through Polya resampling;
  + *Region-wise*: a new pooled sample is created using Polya resampling and then the simulated sites are drawn from this new sample.

All methodological details are found in the article cited above.

See Help.md for a description of the arguments and examples on how to use the function.
