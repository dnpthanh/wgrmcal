# wgrmcal
GRM calculator for MTG2

# Installation
```r
install.packages("devtools")
library(devtools)
devtools::install_github("dnpthanh/wgrmcal")
```

# Example
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
