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
#' data(test_weight)
#' grm_cal(test_freq, test_geno, test_weight)
#' @export
grm_cal <- function(freq_table, geno_table, weight = NULL) {
    # Frequency file
    if (any(is.na(freq_table$MAF))) {
        stop("There is NA(s) in frequency file.")
    }

    freq <- freq_table
    
    # Genotype file
    if (any(is.na(geno_table))) {
        stop("There is NA(s) in genotype file. Imputation may be needed.")
    }

    geno_snp_ids <- colnames(geno_table)[7:ncol(geno_table)]
    common_snp <- intersect(geno_snp_ids, freq$SNP)

    if (length(common_snp) == 0) {
        stop("No common SNPs between freq and geno tables. Please check SNP IDs.")
    }

    freq_matched <- freq[match(common_snp, freq$SNP), ]
    geno_matched <- geno_table[, c(1:6, match(common_snp, colnames(geno_table)))]

    n <- nrow(geno_matched) #number of animals
    p <- dim(geno_matched)[2] - 6 #number of SNPs

    # Prepare M matrix (adjusted -2p)
    geno_adj <- as.matrix(geno_matched[, 7:ncol(geno_matched)]) - 
                2 * matrix(freq_matched$MAF, nrow = n, ncol = p, byrow = TRUE)

    # Weights (default = 1)
    if (is.null(weight)){
        wg <- data.frame(SNP = colnames(geno_matched)[7:dim(geno_matched)[2]], weight = rep(1, p))
    }else{
        w_table <- weight
        cname <- colnames(geno_matched)[7:dim(geno_matched)[2]]
        wg <- data.frame(SNP = cname, weight = rep(1,length(cname)))
        w_ind <- match(w_table[,1], cname)
        valid_idx <- which(!is.na(w_ind))
        wg[w_ind[valid_idx],2] <- w_table[,2][valid_idx]
    }

    # Denominator
    bot <- 2 * sum(freq_matched$MAF * (1 - freq_matched$MAF))

    # MDM'
    grm_matrix <- (as.matrix(geno_adj) %*% diag(wg[,2]) %*% t(as.matrix(geno_adj))) / bot

    # Upper triangle
    idx <- which(upper.tri(grm_matrix, diag = TRUE), arr.ind = TRUE)
    
    # Result
    grm <- data.frame(
        ID_1 = idx[, 2],
        ID_2 = idx[, 1],
        Rela = grm_matrix[idx]
    )

    return(grm)
}
