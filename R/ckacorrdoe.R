#' To Identify Outliers in Design of Experiments with Autocorrelated Errors under lag 1
#'
#' @param trt  Numeric or complex vector containing the treatment levels
#' @param Rep  Numeric or complex vector containing the Replication
#' @param resp Numeric or complex vector containing response variable
#' @param noutli a  number indicating the number of outliers
#'
#' @importFrom stats xtabs acf shapiro.test var cov lm pchisq sd IQR model.matrix
#'
#' @return The output contains Shapiro-Wilk Normality test and Bartlett test for
#' Residuals of the of the model and Cook's distance for each treatment or a
#' combination of treatments in Design of Experiments with auto-correlated errors under lag 1
#'
#' @examples
#' data(ex3)
#' attach(ex3)
#' ckacorrdoe(ex3$trt,ex3$rep,ex3$`sim 1`, 1)
#'
#' @export

ckacorrdoe <- function(trt,Rep,resp,noutli){
  argmts <- lapply(as.list(match.call())[-1], eval)
  missing_vals <- sapply(argmts, function(x) any(is.na(x)))
  ## checking for any missing values
  if(any(missing_vals)){
    warning("One or more arguments contain missing values.")
  }
  cat('------------------------------------------------------------------------\n')
  ## checking for normality assumption
  trt<-as.factor(trt)
  Rep<-factor(Rep)
  model<-lm(resp~Rep+trt)
  n_test<-shapiro.test(model$residuals)
  n_result<-n_test$p.value
  if (n_result > 0.05){
    cat('Normality Assumption : \n')
    print("Normality Assumption is not violated")
  } else {
    cat('Normality Assumption : \n')
    print("Normality Assumption is  violated")
  }
  cat('\n------------------------------------------------------------------------\n\n')
  cat('------------------------------------------------------------------------\n')
  ## check for homogenity of error variances
  bartlett <-function(trt, resp, t, r){
    vari<-tapply(resp,trt,var)
    S2p<-sum((r-1)*vari)/(length(resp)-t)
    A<-(length(resp)-t)*log(S2p)-sum((r-1)*log(vari))
    B<-(1/(3*(t-1)))*(sum(1/(r-1))-(1/(length(resp)-t)))
    Xc1<-A/(1+B)
    pvalue<-1-pchisq(Xc1, t-1)
    output <- pvalue
    return(output)
  }
  barout <- bartlett(trt,resp,length(unique(trt)),length(unique(Rep)))
  if (barout > 0.05){
    cat('Homogenity of variances : \n')
    print("At 5% of significance, residuals can  be considered homocedastic!")
  } else {
    cat('Homogenity of variances : \n')
    print("At 5% of significance, residuals can not be considered homocedastic!")
  }
  cat('\n------------------------------------------------------------------------\n\n')
  ## labeling
  n <- length(trt)
  I_mat <- diag(n)
  r <- length(unique(Rep)) ## number of replications
  v <- length(unique(trt)) ## number of treatments
  j.vec <-  matrix(1, nrow = n, ncol = 1)

  ##generating P matrix
  generate_contrast_matrix <- function(v) {
    if (v < 2) {
      stop("The number of treatments (v) must be at least 2.")
    }

    # Create an empty matrix for (v-1) x v
    P <- matrix(0, nrow = v - 1, ncol = v)

    # Fill the matrix using Helmert contrasts
    for (i in 1:(v - 1)) {
      P[i, 1:i] <- 1
      P[i, i + 1] <- -i
    }

    # Normalize each row (optional)
    P <- t(apply(P, 1, function(row) row / sqrt(sum(row^2))))

    return(P)
  }

  ####
  P <- generate_contrast_matrix(v)

  ##design matrices
  des_trr <- unname(model.matrix(~0+as.factor(trt)))
  des_repr <- unname(model.matrix(~0+as.factor(Rep)))
  one_vecr <- matrix(1, nrow = n, ncol = 1)
  des_fulr <- cbind(one_vecr,des_trr,des_repr)
  omg <- list()
  ## calculation of auto - correlation
  acf_res <- matrix(0, nrow = length(unique(Rep)), ncol = 1)
  for(j in 1:length(unique(Rep))){
    klm <-  xtabs(resp~trt+Rep)
    kat=acf(klm[,j], lag.max = NULL, plot = FALSE)
    acf_res[j] <- kat$acf[2]
  }

  omg <- list()
  for (k in 1:nrow(acf_res)){
    omg[[k]] <- matrix(0, nrow = length(unique(trt)), ncol =length(unique(trt)))##no of trt
  }

  for (k in 1:nrow(acf_res)){
    for (i in 1:length(unique(trt))){
      for (j in 1:length(unique(trt))){
        rho <- (1/(1-acf_res[k]^2))
        if (i==j){
          omg[[k]][i,j]=1
        }
        if (i < j){
          omg[[k]][i,j] = acf_res[k]^(j-i)
        }
        if (i >j) {
          omg[[k]][i,j] = acf_res[k]^(i-j)
        }
      }
    }
    omg[[k]] <- rho * omg[[k]]
  }

  ## diagonally combining all the omg matrices
  omg_fl <- as.matrix(Matrix::bdiag(omg))

  ## calculation of c and q matrices
  bmat_ac <- pracma::inv(omg_fl) - pracma::inv(omg_fl) %*% des_repr %*% MASS::ginv(t(des_repr) %*% pracma::inv(omg_fl) %*% des_repr) %*% t(des_repr) %*% pracma::inv(omg_fl)

  cmat_ac <- t(des_trr) %*% bmat_ac %*% des_trr
  qmat_ac <- t(des_trr) %*% bmat_ac %*% resp

  ## requirements for cook distance
  Vmat_ac <- omg_fl %*% (bmat_ac - bmat_ac %*% des_trr %*% MASS::ginv(cmat_ac) %*% t(des_trr) %*% bmat_ac)
  Hstr <- pracma::inv(omg_fl) %*% Vmat_ac

  ## for cook distance

  ## function for generating trt comb and u_mat
  explore_combinations <- function(data, delete_count) {
    # Get all non-NA indices
    non_na_indices <- which(!is.na(data), arr.ind = TRUE)

    # Generate all combinations of the specified size
    combinations <- utils::combn(seq_len(nrow(non_na_indices)), delete_count)

    # Store the combinations explored
    explored <- list()

    # Iterate through each combination
    for (i in seq_len(ncol(combinations))) {
      # Get the indices for this combination
      combo_indices <- non_na_indices[combinations[, i], , drop = FALSE]

      # Save this combination
      explored[[i]] <- combo_indices
    }
    trt_combs <- data.frame(matrix(nrow = length(explored), ncol = delete_count))
    for (i in seq_along(explored)) {
      for (j in seq_len(nrow(explored[[i]]))) {
        trt_combs[i, j] <- paste0("(", explored[[i]][j, "row"], ", ", explored[[i]][j, "col"], ")")
      }
    }

    # Create matrices marking combinations as 1
    marked_matrices <- list()

    for (i in seq_along(explored)) {
      marked_matrix <- matrix(0, nrow = nrow(data), ncol = ncol(data))
      for (j in seq_len(nrow(explored[[i]]))) {
        marked_matrix[explored[[i]][j, "row"], explored[[i]][j, "col"]] <- 1
      }
      marked_matrices[[i]] <- marked_matrix
    }

    return(list(trt_combs = trt_combs, marked_matrices = marked_matrices))
  }

  ## exp combinations output
  datac_mat <- xtabs(resp~trt+Rep)

  expac_output <- explore_combinations(matrix(datac_mat, nrow = v,ncol = r), delete_count = noutli)
  tr_combac <- expac_output$trt_combs
  markac_mat <- expac_output$marked_matrices

  ## generating u matrix
  u_lis <- list()
  for(i in seq_along(expac_output$marked_matrices)){
    u_vec <- as.vector(markac_mat[[i]])
    u_lis[[i]] <- u_vec
  }

  ckd_acmat <- matrix(0,nrow =length(u_lis), ncol = 1)

  ## loop for cook's distance
  for(outl in 1:length(u_lis)){
    uVu_ac <- t(as.vector(u_lis[[outl]])) %*% pracma::inv(omg_fl) %*% Vmat_ac %*% as.vector(u_lis[[outl]])
    P_thdiff <- (P %*% pracma::pinv(cmat_ac) %*% t(des_trr) %*% bmat_ac %*% as.vector(u_lis[[outl]]) %*% pracma::inv(uVu_ac) %*%
                   t(as.vector(u_lis[[outl]])) %*% pracma::inv(omg_fl) %*% Vmat_ac %*% resp)
    pcp_ac <- P %*% cmat_ac %*% t(P)
    ckd_ac_num <- t(P_thdiff) %*% pcp_ac %*% P_thdiff
    ckd_ac_deno <- var(resp) * (v-1)
    ckd_ac <- ckd_ac_num/ckd_ac_deno
    ckd_acmat[outl,] <- round(ckd_ac,4)
  }

  ## to mark values with threshold values
  ## function to mark
  mark_values_with_asterisk <- function(mat, threshold) {
    # Convert to data frame if not already
    if (!is.data.frame(mat)) {
      mat <- as.data.frame(mat)
    }

    # Rename the column for clarity
    colnames(mat) <- "Cook s Distance"

    # Add a new column to mark values greater than the threshold with an asterisk indicator
    mat$` ` <- ifelse(mat$`Cook s Distance` > threshold, "*", "")

    return(mat)
  }
  fthdac <- IQR(ckd_acmat)
  trshld_ac <- (7/2) * fthdac
  ckd_acmatt <- mark_values_with_asterisk(ckd_acmat,trshld_ac)

  tr_comb2ac <- as.data.frame(tr_combac)
  for (i in seq_along(tr_comb2ac)) {
    colnames(tr_comb2ac)[i] <- paste("(trt_rep)", i, sep = "_")
  }
  ckd_acdf2 <- cbind.data.frame(tr_comb2ac,ckd_acmatt)

  cat('------------------------------------------------------------------------
      \nCook s Distance : \n')
  print(ckd_acdf2)
  cat('\n------------------------------------------------------------------------\n\n')
  cat('------------------------------------------------------------------------
      \nTreatments Suspected as Outliers : \n')
  print(ckd_acdf2[ckd_acdf2[[noutli+2]]== "*",])
}
