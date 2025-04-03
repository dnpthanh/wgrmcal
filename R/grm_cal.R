#' Calculate GRM (Genetic Relationship Matrix)
#'
#' This function computes a genetic relationship matrix based on genotype and frequency data.
#'
#' @param freq_table A data frame containing allele frequencies.
#' @param geno_table A data frame containing genotype data.
#' @param weight_vec A numeric vector of weights (default: NULL).
#' @return A computed GRM matrix.
#' @examples
#' data(test_geno)
#' data(test_freq)
#' grm_cal(test_freq, test_geno)
#' @export
grm_cal <- function(freq_table, geno_table, weight_vec = NULL) {
    # Frequency file
    freq <- freq_table[!is.na(freq_table$MAF), ]
    
    # Genotype file
    geno <- geno_table
    n <- nrow(geno)
    m <- ncol(geno)
    p <- m - 6  # number of SNP

    # Prepare M matrix (adjusted -2p)
    geno_adj <- geno[, 7:m] - 2 * matrix(freq$MAF, nrow = n, ncol = p, byrow = TRUE)

    # Weights (default = 1)
    wg <- if (!is.null(weight_vec)) weight_vec else rep(1, p)

    # Denominator
    bot <- 2 * sum(freq$MAF * (1 - freq$MAF))

    # MDM'
    grm_matrix <- (as.matrix(geno_adj) %*% diag(wg) %*% t(as.matrix(geno_adj))) / bot

    # Upper triangle
    idx <- which(upper.tri(grm_matrix, diag = TRUE), arr.ind = TRUE)
    
    # Result
    grm <- data.frame(
        ID_1 = geno[idx[, 1], 2],
        ID_2 = geno[idx[, 2], 2],
        Rela = grm_matrix[idx]
    )

    return(grm)
}
