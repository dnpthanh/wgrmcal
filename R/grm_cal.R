#' Calculate GRM (Genetic Relationship Matrix)
#'
#' @param freq_table A data frame containing allele frequencies.
#' @param geno_table A data frame containing genotype data.
#' @param weight_vec A numeric vector of weights (default: NULL).
#' @return A computed GRM matrix.
#' @export

grm_cal <- function(freq_table, geno_table, weight_vec = NULL){

    #Frequency file
    freq <- freq_table
    freq2 <- freq[!is.na(freq$MAF),]
    freq2$`1-MAF` <- 1-freq2$MAF

    #Genotype file (raw format)
    geno <- geno_table
    #cname <- colnames(geno)[7:dim(geno)[2]]
    #cname2 <- gsub("_.", "", cname)
    #colnames(geno)[7:dim(geno)[2]] <- cname2

    #Sizes
    n <- nrow(geno)
    m <- ncol(geno)
    p <- m - 6  # number of SNPs
    
    #Prepare for numerator calculation. (Geno - 2p)  
    geno_adj <- geno[, 7:m]  - 2 * matrix(freq2$MAF, nrow = n, ncol = p, byrow = TRUE)

    if (!is.null(weight_vec)){
        wg <- weight_vec
    }else{
        wg <- rep(1, p) 
    }

    #Denominator calculation
    bot <- 2 * sum((freq2$MAF)*(freq2$`1-MAF`))


    #Generate preallocated data frame
    num_pairs <- (n * (n + 1)) / 2      # (i, j) in upper triangular 
    grm <- data.frame(
        ID_1 = character(num_pairs),
        ID_2 = character(num_pairs),
        Rela = numeric(num_pairs),
        stringsAsFactors = FALSE
    )

    #Relationship value calculation
    x <- 1
    for (i in 1:n) {
        for (j in i:n) {
            top <- sum(geno_adj[i, ] * geno_adj[j, ] * wg)   # Vectorize multiplication
            grm$Rela[x] <- top / bot
            grm$ID_1[x] <- geno[i, 2]
            grm$ID_2[x] <- geno[j, 2]
            x <- x + 1  
        }
    }
    return(grm)
}
