#' Calculate GCTA-like standardized GRM (Genetic Relationship Matrix)
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
#' data(test_weight)
#' grm_cal(test_freq, test_geno, test_weight)
#' @export
grm_cal_sd <- function(freq_table, geno_table, weight = NULL) {
    # Frequency file
    if (any(is.na(freq_table$MAF))) {
        stop("There is NA(s) in frequency file.")
    }
    
    # Genotype file
    if (any(is.na(geno_table))) {
        stop("There is NA(s) in genotype file. Imputation may be needed.")
    }

    freq <- freq_table
    geno_snp_ids <- colnames(geno_table)[7:ncol(geno_table)]
    common_snp <- intersect(geno_snp_ids, freq$SNP)

    if (length(common_snp) == 0) {
        stop("No common SNPs between frequency and genotype tables. Please check SNP IDs.")
    }

    freq_matched <- freq[match(common_snp, freq$SNP), ]
    geno_matched <- geno_table[, c(1:6, match(common_snp, colnames(geno_table)))]

    n <- nrow(geno_matched) # number of animals
    p <- ncol(geno_matched) - 6 # number of SNPs

    # Genotype matrix
    geno_raw <- as.matrix(geno_matched[, 7:ncol(geno_matched)])
    maf <- freq_matched$MAF

    # Standardize: (G - 2p) / sqrt(2p(1-p))
    geno_adj <- sweep(geno_raw, 2, 2 * freq_matched$MAF, "-") / sqrt(2 * maf * (1 - maf))

    # Weight
    if (is.null(weight)) {
        wg <- rep(1, p)
    } else {
        cname <- colnames(geno_matched)[7:ncol(geno_matched)]
        wg <- rep(1, length(cname))
        w_ind <- match(weight[, 1], cname)
        valid_idx <- which(!is.na(w_ind))
        wg[w_ind[valid_idx]] <- weight[, 2][valid_idx]
    }

    # GRM = ZZ' / m
    ZW <- geno_adj %*% diag(wg)
    grm_matrix <- tcrossprod(ZW) / sum(wg)  # if having weights, divide by sum of weights

    # Upper tri
    idx <- which(upper.tri(grm_matrix, diag = TRUE), arr.ind = TRUE)

    # Output
    grm <- data.frame(
        ID_1 = idx[, 2],
        ID_2 = idx[, 1],
        Rela = grm_matrix[idx]
    )

    return(grm)
}
