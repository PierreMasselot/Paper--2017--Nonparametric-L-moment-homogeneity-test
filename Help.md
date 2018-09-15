
# Lmom.homog.test
## Nonparametric L-moment honogeneity test}

### Description
  Performs the Nonparametric Hosking-Wallis L-moment homogeneity test. Includes the parametric procedure of Hosking and Wallis and several nonparametric procedures.  
  
### Usage
       Lmom.homog.test(X, Nsim = 500, type = "B", typenorm = "2", center = type %in% c("B","Yr","Ys"), printInfo = T)

### Arguments 
  **X**: List. Each element of the list contains the data for one site. Each element is a matrix with one variable per column (including the univariate case). This allows different record length for each site
  
  **Nsim**: Numeric. Number of simulated regions for the estimation of the V distribution
  
  **type**: Character. type of simulation to perform (HW : classical Kappa, M : permutations, B : bootstrap, Ys : Polya site-wise, Yr : polya region-wise)
  
  **typenorm**: Type of norm to use (see the function [norm](http://127.0.0.1:20740/library/base/html/norm.html))
  
  **center**: Logical. Center the sites? Recommended for the Bootstrap test
  
  **printInfo**: Logical. If TRUE (the default), informations about the test and simulations iterations are displayed in the console

### Details
  The test is based on L-moments, which are linear combinations of order statistics and are analogous to classical moments. It is intended to test the homogeneity between a set of sites, up to a scale factor. Thus, the test statistic `V` is the weighted mean of the second order L-moments of each sample.
  
  The distribution of `V` under the null hypothesis of homogeneity is estimated through simulation. The type of simulation is determined by the argument `type`. When `type = 'HW'`, a four parameter Kappa distribution is fit to the pooled sample and simulations are drawn from the fitted distribution. All other values of `r type` correspond to non parametric simulation. `type = 'B'` (the default) corresponds to bootstrap, `type = 'M'` corresponds to permutations, `type = 'Ys'` to polya resampling used to directly draw the sites and `type = 'Yr'` to Polya resampling used *a priori* to simulate a new sample from which sites are drawn. 

### Value
  A list containing the components:
  
  - **Vobs**: The observed V statistic
  
  - **Vsim**: A vector of length Nsim containing all simulated V statistics
  
  - **H**: The final H statistic
  
  - **p.value**: The p-value of the test 

### References
  Hosking, J., Wallis, J., 1993. Some statistics useful in regional frequency analysis. *Water Resources Research* 29, 271–281.
  
  Masselot, P., Chebana, F., Ouarda, T.B.M.J., 2017. Fast and direct nonparametric procedures in the L-moment homogeneity test. *Stochastic Environmental Research and Risk Assessment* 31, 509–522.

### Examples
  ```
  library(lmomRFA)
  data(Maxwind)

  # Bootstrap simulations
  resB <- Lmom.homog.test(Maxwind)

  # Permutations
  resM <- Lmom.homog.test(Maxwind, type = "M")
  ```
