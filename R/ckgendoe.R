#' To Identify Outliers in Design of Experiments
#'
#' @param trt Numeric or complex vector containing the treatment levels
#' @param Rep Numeric or complex vector containing the Replication
#' @param resp Numeric or complex vector containing response variable
#' @param noutli a number indicating the number of outliers
#'
#' @importFrom stats xtabs acf shapiro.test var cov lm pchisq sd IQR model.matrix
#'
#' @return The output contains Shapiro-Wilk Normality test and Bartlett test for
#' Residuals of the of the model and Cook's distance for each treatment or a
#' combination of treatments
#'
#' @export
#'
#' @examples
#' # example code
#' data(ex1)
#' ckgendoe(ex1$trt,ex1$rep,ex1$yld,1)


ckgendoe <- function(trt,Rep,resp,noutli){
  argmts <- lapply(as.list(match.call())[-1], eval)
  missing_vals <- sapply(argmts, function(x) any(is.na(x)))
  ## checking for any missing values
  if(any(missing_vals)){
    warning("One or more arguments contain missing values.")
  }
  cat('------------------------------------------------------------------------\n')
  ## checking for normality assumption
  trt<-as.factor(trt)
  Rep<-as.factor(Rep)
  model<-lm(resp~Rep+trt)
  n_test<-shapiro.test(model$residuals)
  n_result<-n_test$p.value
  if (n_result > 0.05){
    cat('Normality Assumption : \n')
    print("At 5% of significance, Normality Assumption is not violated")
  } else {
    cat('Normality Assumption : \n')
    print("At 5% of significance, Normality Assumption is  violated")
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

  ##design matrices
  des_trr <- unname(model.matrix(~0+as.factor(trt)))
  des_repr <- unname(model.matrix(~0+as.factor(Rep)))
  one_vecr <- matrix(1, nrow = n, ncol = 1)
  des_fulr <- cbind(one_vecr,des_trr,des_repr)

  ## calculation of c and q matrices
  cmat <- (t(des_trr) %*% des_trr) - ((t(des_trr) %*% des_repr) %*% MASS::ginv(t(des_repr) %*% des_repr) %*%
                                        (t(des_repr) %*% des_trr))
  qmat <- (t(des_trr) %*% resp) - ((t(des_trr) %*% des_repr) %*% MASS::ginv(t(des_repr) %*% des_repr) %*%
                                     (t(des_repr) %*% resp))
  bmat <- I_mat - des_repr %*% MASS::ginv(t(des_repr) %*% des_repr) %*% t(des_repr)

  ## function to obtain trt combinations and U matrix
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

    # Create trt_combs
    trt_combs <- data.frame(matrix(nrow = length(explored), ncol = delete_count))
    for (i in seq_along(explored)) {
      for (j in seq_len(nrow(explored[[i]]))) {
        trt_combs[i, j] <- paste0("(", explored[[i]][j, "row"], ", ", explored[[i]][j, "col"], ")")
      }
    }

    # Create nested marked_matrices structure
    marked_matrices <- vector("list", length = length(explored))

    for (i in seq_along(explored)) {
      # Initialize a list for this combination
      combination_list <- vector("list", length = delete_count)

      for (j in seq_len(nrow(explored[[i]]))) {
        # Create a matrix initialized with zeros
        marked_matrix <- matrix(0, nrow = nrow(data), ncol = ncol(data))
        # Mark the specific treatment as 1
        marked_matrix[explored[[i]][j, "row"], explored[[i]][j, "col"]] <- 1
        # Save the matrix in the sub-list for this combination
        combination_list[[j]] <- marked_matrix
      }

      # Save the sub-list under the corresponding combination
      marked_matrices[[i]] <- combination_list
    }

    # Output: List of all trt_combs and nested marked_matrices
    return(list(trt_combs = trt_combs, marked_matrices = marked_matrices))
  }

  ##explore combinations function output
  data_fmat <- xtabs(resp~trt+Rep)
  exp_out <- explore_combinations(matrix(data_fmat, nrow = v,ncol = r), delete_count = noutli)

  ## trt combinations
  tr_comb <- exp_out$trt_combs
  mark_mats <- exp_out$marked_matrices

  ## now generating u matrix for given t
  U_mat <- list()
  for(i in seq_along(mark_mats)){
    u_ini_mat <- matrix(0,nrow = n, ncol = noutli)
    for(j in 1:noutli){
      u_vec <- as.vector(mark_mats[[i]][[j]])
      u_ini_mat[,j] <- u_vec
    }
    U_mat[[i]] <- u_ini_mat
  }

  ## requirements for cook distance
  S_mat <- bmat %*% des_trr %*% pracma::pinv(cmat) %*% t(des_trr) %*% bmat
  V_mat <- bmat - S_mat

  ## for cook distance
  ckd_mat <- matrix(0,nrow = length(U_mat), ncol = 1)

  for (i in 1:length(U_mat)){
    hdel <- MASS::ginv(t(U_mat[[i]]) %*% V_mat %*% U_mat[[1]]) %*% t(U_mat[[i]]) %*% V_mat %*% resp
    cknum <- t(hdel) %*% t(U_mat[[i]]) %*% S_mat %*% U_mat[[i]] %*% hdel
    ckden <- (v-1) * var(resp)
    ckdis <- cknum/ckden
    ckd_mat[i,] <-  round(ckdis,4)
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
  fthd <- IQR(ckd_mat)
  trshld <- (7/2) * fthd
  ckd_matt <- mark_values_with_asterisk(ckd_mat,trshld)

  tr_comb2 <- as.data.frame(tr_comb)
  for (i in seq_along(tr_comb2)) {
    colnames(tr_comb2)[i] <- paste("(trt_rep)", i, sep = "_")
  }
  ckd_df2 <- cbind.data.frame(tr_comb2,ckd_matt)

  cat('------------------------------------------------------------------------
      \nCook s Distance : \n')
  print(ckd_df2)
  cat('\n------------------------------------------------------------------------\n\n')
  cat('------------------------------------------------------------------------
      \nTreatments Suspected as Outliers : \n')
  print(ckd_df2[ckd_df2[[noutli+2]]== "*",])
}
