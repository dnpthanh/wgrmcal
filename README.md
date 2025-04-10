# wgrmcal
GRM calculator for MTG2

# Requirements
- Genotype file must not have NA(s).
- Imputation may be used for imputing missing genotype.
- SNP IDs must be matched between genotype file and frequency file.
- Parallel is recommended before executing grm_cam() function.

# Installation
```r
install.packages("devtools")
library(devtools)
devtools::install_github("dnpthanh/wgrmcal")
```

# Example

```r
# Optional

library(parallel)
numCores <- detectCores()  # number of CPU cores
cl <- makeCluster(numCores - 1)  # Use maximum number of CPU cores (changable)

library(doParallel)
registerDoParallel(cores = numCores - 1)
```

```r
library(wgrmcal)
```

```r
# without weights
data(test_geno)
data(test_freq)
test_grm <- grm_cal(test_freq, test_geno)
```

```r
# with weights
w <- rep(0.5, 14)
test_grm2 <- grm_cal(test_freq, test_geno, w)
```
