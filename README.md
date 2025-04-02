# wgrmcal
GRM calculator for MTG2

# example

# without weights
test_geno <- data(test_geno)
test_freq <- data(test_freq)
test_grm <- grm_cal(test_freq, test_geno)

# with weights
w <- rep(0.5, 14)
test_grm2 <- grm_cal(test_freq, test_geno, w)
