#' @import coloc 
#' @import data.table

#'@param t1 SuSiEx .snp file output from trait 1
#'@param t2 SuSiEx .snp file output from trait 2
#'@param p1 prior probability that a SNP is causal for trait 1
#'@param p2 prior probability that a SNP is causal for trait 2
#'@param p12 prior probability that a SNP is causal for both traits
#'@return coloc style result for coloc_SuSiEx
#'@export


coloc_susiex <- function(t1, t2, p1 = 1e-4, p2=1e-4, p12=5e-6, overlap.min=0.5, trim_by_posterior=T) {
  
  t1 <- as.data.frame(t1)
  t2 <- as.data.frame(t2)
  
  ncst1 <- unique(gsub(".*\\((CS[0-9]+).*", "\\1", names(t1)[grep("LogBF\\(CS", names(t1))]))
  t1matrix <- matrix(ncol = nrow(t1), nrow = length(ncst1))
  
  ncst2 <- unique(gsub(".*\\((CS[0-9]+).*", "\\1", names(t2)[grep("LogBF\\(CS", names(t2))]))
  t2matrix <- matrix(ncol = nrow(t2), nrow = length(ncst2))
  
  if ((length(ncst1) ==0) | (length(ncst2) == 0)) {
    return(NA)
  }
  
  for (k in 1:length(ncst1)) {
    cs_col <- grep(paste0(ncst1[k], ",Pop"), names(t1), value = TRUE)
    sumlogbf <- rowSums(t1[,which(colnames(t1)%in% cs_col)])
    t1matrix[k,] <- sumlogbf
  }
  
  colnames(t1matrix) <- t1$SNP
  
  for (k in 1:length(ncst2)) {
    cs_col <- grep(paste0(ncst2[k], ",Pop"), names(t2), value = TRUE)
    sumlogbf <- rowSums(t2[,which(colnames(t2)%in% cs_col)])
    t2matrix[k,] <- sumlogbf
  }
  
  colnames(t2matrix) <- t2$SNP
  
  coloc_susiex_result <- coloc:::coloc.bf_bf(t1matrix, t2matrix, p1, p2, p12, overlap.min, trim_by_posterior)
  
  return(coloc_susiex_result)
  
}

#'@param t1 mscaviar post.txt output file for trait 1
#'@param t2 mscaviar post.txt output file for trait 2
#'@return msecaviar colocalization results
#'@export

msecaviar <- function(t1,t2) {
  clpp <- t1$Prob_in_pCausalSet * t2$Prob_in_pCausalSet
  scaled_clpp <- clpp / sum(clpp)
  locus_level_clpp <- sum(clpp)
  snps <- t1$SNP_ID
  
  msecaviar_result <- list(data.frame(snps, scaled_clpp, clpp), locus_level_clpp)
  names(msecaviar_result) <- c("variant_level", "locus_level_clpp")
  
  return(msecaviar_result)
}

#'@param t1 SuSiEx .snp file output from trait 1
#'@param t2 SuSiEx .snp file output from trait 2
#'@return SuSiEx_eCAVIAR colocalization results
#'@export

susiex_ecaviar <- function(t1, t2) {
  
  t1 <- as.data.frame(t1)
  t2 <- as.data.frame(t2)
  
  pip_cols1 <- grep("PIP", names(t1), value = TRUE)
  pip_cols2 <- grep("PIP", names(t2), value = TRUE)
  
  if ((length(pip_cols1) == 0) | (length(pip_cols2) == 0)) { return(NA)}
  
  t1pip <- as.matrix(t1[, which(colnames(t1) %in% pip_cols1)])
  colnames(t1pip) <- pip_cols1
  t2pip <- as.matrix(t2[, which(colnames(t2) %in% pip_cols2)])
  colnames(t2pip) <- pip_cols2
  
  # t11pip <- t11[,..pip_cols1]
  # t22pip <- t22[,..pip_cols2]
  
  coloc_comparisons <- expand.grid(colnames(t1pip), colnames(t2pip), stringsAsFactors = FALSE)
  # clpp_list <- mapply(function(col1, col2) {
  #   t11pip[[col1]] * t22pip[[col2]]
  #  }, coloc_comparisons[,1], coloc_comparisons[,2], SIMPLIFY = FALSE)

  
  clpp_list <- mapply(function(col1, col2) {
    t1pip[,which(colnames(t1pip) == col1)] * t2pip[,which(colnames(t2pip) == col2)]
  }, coloc_comparisons[,1], coloc_comparisons[,2], SIMPLIFY = FALSE)
  
  names(clpp_list) <- paste(sub(".*\\((.*)\\).*", "\\1", coloc_comparisons[,1]), 
                            "x", sub(".*\\((.*)\\).*", "\\1", coloc_comparisons[,2]), sep = "_")
  
  clpp <- (as.data.frame(clpp_list))
  scaled_clpp <- sweep(clpp, 2, colSums(clpp), FUN = "/")
  locus_level_clpp <- colSums(clpp)
  
  clpp$snps <- t1$SNP
  scaled_clpp$snps <- t1$SNP
  
  
  results <- list(clpp, scaled_clpp, locus_level_clpp)
  names(results) <- c("clpp", "scaled_clpp", "locus_level_clpp")
  
  
  return(results)
  
}

#'@param t1 MsCAVIAR post.txt output from trait 1
#'@param t2 MsCAVIAR post.txt output from trait 2
#'@param p1 prior probability that a SNP is causal for trait 1
#'@param p2 prior probability that a SNP is causal for trait 2
#'@param p12 prior probability that a SNP is causal for both traits
#'@return coloc style result for MsCAVIAR_coloc
#'@export

mscaviar_coloc <- function(t1, t2, p1 = 1e-4, p2=1e-4, p12=5e-6, overlap.min=0.5, trim_by_posterior=T) {
  
  t1$Prob_in_pCausalSet[which(t1$Prob_in_pCausalSet == 1)] <- 1 - 1e-12
  t1$Prob_in_pCausalSet[which(t1$Prob_in_pCausalSet == 0)] <- 1e-12
  t2$Prob_in_pCausalSet[which(t2$Prob_in_pCausalSet == 1)] <- 1 - 1e-12
  t2$Prob_in_pCausalSet[which(t2$Prob_in_pCausalSet == 0)] <- 1e-12
  
  bft1 <- t1$Prob_in_pCausalSet/(1-t1$Prob_in_pCausalSet) * (1-p1)/p1
  logbft1 <- as.vector(unlist(log(bft1)))
  bft2 <- t2$Prob_in_pCausalSet/(1-t2$Prob_in_pCausalSet) * (1-p2)/p2
  logbft2 <- as.vector(unlist(log(bft2)))
  names(logbft1) <- names(logbft2) <- as.vector(unlist(t1$SNP_ID))
  
  mscaviar_coloc_result <- coloc:::coloc.bf_bf(logbft1, logbft2, p1, p2, p12, overlap.min, trim_by_posterior)
  return(mscaviar_coloc_result)
  
}


#'@param cs_res coloc_SuSiEx result 
#'@param rowind coloc comparison of interest to generate a credible set for (which row)
#'@export

make_coloc_susiex_credible_set <- function(cs_res, rowind) {
  cs_res$results <- as.data.frame(cs_res$results)
  cs_res_sorted <- cs_res$results[order(cs_res$results[,rowind], decreasing = TRUE),]
  cumsum <- cumsum(cs_res_sorted[,rowind])
  if (cs_res_sorted[1,rowind] > 0.95) {
    cs <- cs_res_sorted[1,]
  } else {
    cs <- cs_res_sorted[which(cumsum <= 0.95),]
    while (cumsum(cs[,rowind])[nrow(cs)] < 0.95) {
      cs <- rbind(cs, cs_res_sorted[nrow(cs)+1,])
    }
  }
  return(cs)
}


#'@param mse_res msecaviar results
#'@export
make_msecaviar_credible_set <- function(mse_res) {
  mse_res_sorted <- mse_res[order(-mse_res$scaled_clpp),]
  cumsum <- cumsum(mse_res_sorted$scaled_clpp)
  if (cumsum[1] > 0.95) {
    cs <- mse_res_sorted[1,]
  } else {
    cs <- mse_res_sorted[cumsum <= 0.95,]
    if (cumsum(cs[nrow(cs),2]) < 0.95) {
      cs <- rbind(cs, mse_res_sorted[nrow(cs)+1,])
    }
  }
  
  return(cs)
}


#'@param se_res susiex_ecaivar results
#'@param rowind coloc comparison to generate a credible set for (which row)
#'@export
make_susiex_ecaviar_credible_set <- function(se_res, rowind) {
  se_res_sorted <- se_res[order(se_res[,rowind], decreasing = TRUE),]
  cumsum <- cumsum(se_res_sorted[,rowind]) 
  if (se_res_sorted[1,rowind] > 0.95) {
    cs <- se_res_sorted[1,]
  } else {
    cs <- se_res_sorted[which(cumsum <= 0.95),]
    while (cumsum(cs[,rowind])[nrow(cs)] < 0.95) {
      cs <- rbind(cs, se_res_sorted[nrow(cs) + 1,])
    }
  }
  return(cs)
}


#'@param mc_res mscaviar_coloc results
#'@export
make_mscaviar_coloc_credible_set <- function(mc_res) {
  mc_res_sorted <- mc_res$results[order(mc_res$results$SNP.PP.H4.abf, decreasing = TRUE),]
  cumsum <- cumsum(mc_res_sorted$SNP.PP.H4.abf)
  if (mc_res_sorted$SNP.PP.H4.abf[1] > 0.95) {
    cs <- mc_res_sorted[1,]
  } else {
    cs <- mc_res_sorted[cumsum <= 0.95,]
    if (cumsum(cs[nrow(cs),2]) < 0.95) {
      cs <- rbind(cs, mc_res_sorted[nrow(cs)+1,])
    }
  }
}
