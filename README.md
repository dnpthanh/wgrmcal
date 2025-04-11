# wgrmcal
GRM calculator for MTG2
grm_cal() function is used for calculating GRM without standardizing.
grm_cal_sd() function is used for calculating standardized GRM (GCTA-like).

# Requirements
- Genotype file must not have NA(s).
- Imputation should be performed for imputing missing genotype.
- SNP IDs must be matched between genotype file and frequency file.
- Output GRM is in .txt format for MTG2, first and second columns are ID order of fam files.

# Installation
```r
install.packages("devtools")
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
test_grm_sd <- grm_cal_sd(test_freq, test_geno)
```

```r
# with weights
data(test_weight)
test_grm2 <- grm_cal(test_freq, test_geno, test_weight)
test_grm_sd2 <- grm_cal_sd(test_freq, test_geno, test_weight)
```
