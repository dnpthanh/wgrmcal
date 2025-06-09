#' Calculate GRM (Genetic Relationship Matrix)
#'
#' This function computes a genetic relationship matrix based on genotype and frequency data.
#'
#' @param freq_table A data frame containing allele frequencies.
#' @param geno_table A data frame containing genotype data.
#' @param weight A data frame includes SNP names in the first column and their weights in the second column (default: NULL).
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
        stop("There is NA(s) in genotype file. Phenotype may appear as NAs or Imputation may be needed.")
    }

    geno_snp_ids <- colnames(geno_table)[7:ncol(geno_table)]
    common_snp <- intersect(geno_snp_ids, freq$SNP)

    if (length(common_snp) == 0) {
        stop("No common SNPs between frequency and genotype tables. Please check SNP IDs.")
    }

    freq_matched <- freq[match(common_snp, freq$SNP), ]
    geno_matched <- geno_table[, c(1:6, match(common_snp, colnames(geno_table)))]

    n <- nrow(geno_matched) #number of animals
    p <- dim(geno_matched)[2] - 6 #number of SNPs

    # Weights (default = 1)
    snp_names <- colnames(geno_matched)[7:dim(geno_matched)[2]]
    if (is.null(weight)){
        wg <- data.frame(SNP = snp_names, weight = rep(1, p))
    }else{
        wg <- data.frame(SNP = snp_names, weight = rep(1,p))
        matched_idx <- match(weight[, 1], snp_names)
        valid_idx <- which(!is.na(matched_idx))
        wg$weight[matched_idx[valid_idx]] <- weight[, 2][valid_idx]
    }

    # Prepare M matrix (x - 2p)
    geno_adj <- as.matrix(geno_matched[, 7:ncol(geno_matched)]) - 
                2 * matrix(freq_matched$MAF, nrow = n, ncol = p, byrow = TRUE)

    # 2p(1-p) for each SNP
    dk <- 2 * freq_matched$MAF * (1 - freq_matched$MAF)
    
    # Divide each column of geno_adj by 2p(1 - p)
    geno_scaled <- sweep(geno_adj, 2, dk, FUN = "/")
    
    # Calculate GRM
    grm_matrix <- (as.matrix(geno_scaled) %*% diag(wg$weight) %*% t(as.matrix(geno_scaled))) / sum(wg$weight)

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
