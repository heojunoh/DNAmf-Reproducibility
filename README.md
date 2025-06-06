# Diffusion Non-Additive Model for Multi-Fidelity Computer Experiments with Tuning Parameters (Reproducibility)

Junoh Heo, Romain Boutelet, and Chih-Li Sung Jun 6, 2025

This instruction aims to reproduce the results in the paper "*Diffusion Non-Additive Model for Multi-Fidelity Computer Experiments with Tuning Parameters*".

The following results are reproduced in this file

-   Section 4: Figure 3 and 4
-   Section 5: Figure 6 and S1
-   Section 6: Figure 7 and 8

The approximate running times for each section are as follows:

-   Section 4: \~50 minutes
-   Section 5: \~48 hours
-   Section 6: \~2 minutes

$$WARNING$$ Reproducing the results of Section 5 can take a considerable amount of time (more than a day). We recommend reproducing the results of Section 4 and Section 6 first before attempting Section 5.

The package `MuFiCokriging` has been removed from the CRAN repository. You can download `MuFiCokriging` using the following code:

``` r
library(devtools)
install_github("cran/MuFiCokriging")
```

##### Step 0.1: setting

As mentioned earlier, reproducing all the results can be time-consuming. Therefore, we have included options in the settings to allow users to skip the most resource-intensive parts. Specifically, the logical object `FEM` determines whether you want to run the FEM for comparison, which is implemented in Matlab. This involves finite element simulations and can take a long time to run. Therefore, we recommend setting **FEM=FALSE** to avoid running this method during your initial attempt.

``` r
setwd("/DNAmf-Reproducibility/") # set your working directory
eps <- sqrt(.Machine$double.eps) # small nugget for numeric stability
FEM <- TRUE # do you want to run FEM
```

##### Step 0.2: load functions and packages

``` r
install.packages("MuFiMeshGP.AL_1.0.tar.gz", repos = NULL, type = "source") # additional package for comparison
library(lhs)
library(plgp)
library(MuFiCokriging)
library(MuFiMeshGP.AL)
library(RNAmf)
library(matlabr)
library(ggplot2)
library(ggpubr)
library(hrbrthemes)

source("GP.R")
source("GP_nonsep.R")
source("DNAmf_sqex.R")
source("predict.DNAmf_nonsep.R")
source("NestedX.R") 
source("crps.R") # CRPS evaluation
```

## Section 4:

##### Reproducing Figure 3

This reproduces the illustrations of the additive function (top row) and the non-additive function (bottom row).

``` r
# Run code
source("Numerical Studies Figure 3.R")
figure3
```

<img src="figure/Figure 3.png" style="display: block; margin: auto;"/>

##### Reproducing Figure 4

This reproduces the boxplot results for Section 4.

``` r
source("Numerical Studies Additive comparison.R")
source("Numerical Studies Non-additive comparison.R")
source("Numerical Studies Currin comparison.R")
source("Numerical Studies Borehole comparison.R")
source("Numerical Studies Figure 4.R")
figure4
```

<img src="figure/Figure 4.png" style="display: block; margin: auto;"/>

## Section 5

This section includes MATLAB code for running the finite element simulations for solving the Poisson's equation and the vibration of square plate. Before running .m files, you may need to set a proper working environment.

##### Reproducing Figure 6

This reproduces the boxplot results for Section 5. Running FEM can take a significant amount of time. We recommend setting **FEM=FALSE** during your initial attempt to avoid long running times.

``` r
if(FEM) source("Real Applications Poisson comparison.R")
if(FEM) source("Real Applications Plate comparison.R")
source("Real Applications time-dependent PDE comparison.R")
source("Real Applications Figure 6.R")
figure6
```

<img src="figure/Figure 6.png" style="display: block; margin: auto;"/>

##### Reproducing Figure S1

This reproduces the illustrations of Poisson's equation (top row) and the heat equation (bottom row).

``` r
# Run code
source("Real Applications Figure S1.R")
figureS1
```

<img src="figure/Figure S1.png" style="display: block; margin: auto;"/>

## Section 6

To generate pseudo outputs, imputer.R needs to be loaded.

``` r
source("imputer.R")
```

##### Reproducing Figure 7

This produces the illustrations of pseudo-data generation.

``` r
source("DNA Model with Non-nested Design Figure 7.R")
figure7
```

<img src="figure/Figure 7.png" style="display: block; margin: auto;"/>

##### Reproducing Figure 8

This produces the illustrations of the DNA model with non-nested design.

``` r
source("DNA Model with Non-nested Design Figure 8.R")
figure8
```

<img src="figure/Figure 8.png" style="display: block; margin: auto;"/>
