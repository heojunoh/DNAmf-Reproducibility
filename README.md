# Diffusion Non-Additive Model for Multi-Fidelity Computer Experiments with Tuning Parameters (Reproducibility)

Junoh Heo, Romain Boutelet, Wenjia Wang, and Chih-Li Sung Apr 20, 2026

This instruction aims to reproduce the results in the paper "*Diffusion Non-Additive Model for Multi-Fidelity Computer Experiments with Tuning Parameters*".

The following results are reproduced in this file

-   Section 4: Figure 3
-   Section 5: Figure 4, 5, S1, S2, S3, and S4
-   Section 6: Figure 7 and S5
-   Section 7: Figure 8

The approximate running times for each section are as follows:

-   Section 4: \~5 seconds
-   Section 5: \~2.5 hours
-   Section 6: \~30 minutes
-   Section 7: \~2 minutes

The package `MuFiCokriging` has been removed from the CRAN repository. You can download `MuFiCokriging` using the following code:

``` r
library(devtools)
install_github("cran/MuFiCokriging")
```

##### Step 0.1: setting

Reproducing all the results can be time-consuming. Therefore, we have included options in the settings to allow users to skip the most resource-intensive parts. Specifically, the logical object NARGP determines whether you want to run the NARGP method for comparison, which is implemented in Python. This method involves sophisticated Monte Carlo (MC) approximations and can take a long time to run. Additionally, installing the required Python packages can be complex and time-consuming. Therefore, we recommend setting NARGP=FALSE to avoid running this method during your initial attempt.

``` r
setwd("/DNAmf-Reproducibility/") # set your working directory
eps <- sqrt(.Machine$double.eps) # small nugget for numerical stability
NARGP <- TRUE # do you want to run NARGP
```

##### Step 0.2: load functions and packages

``` r
library(lhs)
library(plgp)
library(MuFiCokriging)
library(MuFiMeshGP)
library(RNAmf)
library(matlabr)
library(ggplot2)
library(ggpubr)
library(ReacTran)
library(deSolve)

source("GP.R")
source("GP_nonsep.R")
source("DNAmf_sqex.R")
source("predict.DNAmf_nonsep.R")
source("NestedX.R") 
source("crps.R") # CRPS evaluation

# You may need to install packages using pip install on terminal or command. Sometimes it requires to restart R/Rstudio after packages installed.
if(NARGP){
  library(reticulate)
  py_install("GPy==1.10.0")
  py_install("numpy==1.24.0")
  py_install("pandas==2.2.1")
  py_install("scipy==1.13.0")
  py_install("rpy2==3.5.11")
}

# Moreover, MATLAB should be appropriately installed to run the results(Poisson and Plate) in Section 6.
```

## Section 4

##### Reproducing Figure 3

This produces the illustrations of the optimal multi-fidelity design.

``` r
source("Optimal multi-fidelity design Figure 3.R")
figure3
```

<img src="figure/Figure 3.png" style="display: block; margin: auto;"/>

## Section 5:

##### Reproducing Figure 4

This reproduces the illustrations of the additive function (top row) and the non-additive function (bottom row).

``` r
# Figure 4
source("Numerical Studies Figure 4.R")
figure_illustration
```

<img src="figure/Figure 4.png" style="display: block; margin: auto;"/>

##### Reproducing Figure 5, S1, S2, S3, and S4

This reproduces the boxplot results for Section 5 Numerical Studies.

``` r
source("GP_nonsep.R")
# Figure 5
c <- 0.7; gam <- -2
source("Numerical Studies Additive comparison.R")
source("Numerical Studies Non-additive comparison.R")
source("Numerical Studies Currin comparison.R")
source("Numerical Studies Borehole comparison.R")
source("Numerical Studies Figure 5.R")
figure_numerical
```

<img src="figure/Figure 5.png" style="display: block; margin: auto;"/>

This reproduces the boxplot results for Supplementary Materials with c = 0.75 and gamma = 2.

``` r
# Figure S1
c <- 0.75; gam <- -2
source("Numerical Studies Additive comparison.R")
source("Numerical Studies Non-additive comparison.R")
source("Numerical Studies Currin comparison.R")
source("Numerical Studies Borehole comparison.R")
source("Numerical Studies Figure 5.R")
figure_numerical
```

<img src="figure/Figure S1.png" style="display: block; margin: auto;"/>

This reproduces the boxplot results for Supplementary Materials with c = 0.5 and gamma = 1.

``` r
# Figure S2
c <- 2/3; gam <- -1.5
source("Numerical Studies Additive comparison.R")
source("Numerical Studies Non-additive comparison.R")
source("Numerical Studies Currin comparison.R")
source("Numerical Studies Borehole comparison.R")
source("Numerical Studies Figure 5.R")
figure_numerical
```

<img src="figure/Figure S2.png" style="display: block; margin: auto;"/>

This reproduces the boxplot results for Supplementary Materials with c = 2/3 and gamma = 1.5.

``` r
# Figure S3
c <- 0.55; gam <- -1
source("Numerical Studies Additive comparison.R")
source("Numerical Studies Non-additive comparison.R")
source("Numerical Studies Currin comparison.R")
source("Numerical Studies Borehole comparison.R")
source("Numerical Studies Figure 5.R")
figure_numerical
```

<img src="figure/Figure S3.png" style="display: block; margin: auto;"/>

This reproduces the boxplot results for Supplementary Materials with c = 0.55 and gamma = 1.

``` r
# Figure S4
c <- 0.5; gam <- -1
source("Numerical Studies Additive comparison.R")
source("Numerical Studies Non-additive comparison.R")
source("Numerical Studies Currin comparison.R")
source("Numerical Studies Borehole comparison.R")
source("Numerical Studies Figure 5.R")
figure_numerical
```

<img src="figure/Figure S4.png" style="display: block; margin: auto;"/>

## Section 6

##### Reproducing Figure 7

This reproduces the boxplot results for Section 6.

``` r
source("Real Applications Poisson comparison.R")
source("Real Applications Plate comparison.R")
source("Real Applications time-dependent PDE comparison.R")
source("Real Applications Figure 7.R")
figure7
```

<img src="figure/Figure 7.png" style="display: block; margin: auto;"/>

##### Reproducing Figure S5

This reproduces the illustrations of Poisson's equation (top row) and the heat equation (bottom row).

``` r
# Run code
source("Real Applications Figure S5.R")
figureS5
```

<img src="figure/Figure S5.png" style="display: block; margin: auto;"/>

## Section 7

To generate pseudo outputs, imputer.R needs to be loaded.

``` r
source("imputer.R")
```

##### Reproducing Figure 8

This produces the illustrations of the DNA model with non-nested design.

``` r
source("DNA Model with Non-nested Design Figure 8.R")
figure8
```

<img src="figure/Figure 8.png" style="display: block; margin: auto;"/>
