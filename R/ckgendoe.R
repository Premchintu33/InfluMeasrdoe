#' To Identify Outliers in Design of Experiments
#'
#' @param trt Numeric or complex vector containing the treatment levels
#' @param Rep Numeric or complex vector containing the Replication
#' @param resp Numeric or complex vector containing response variable
#' @param noutli a number indicating the number of outliers
#'
#' @return The output contains Shapiro-Wilk Normality test and Bartlett test for
#' Residuals of the of the model and Cook's distance for each treatment or a
#' combination of treatments
#'
#'
#' @examples
ckgdoe <- function(trt,Rep,resp,noutli){
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
  v <- length(unique(trt)) ## number of treatments
  j.vec <-  matrix(1, nrow = n, ncol = 1)

  ##design matrices
  des_trr <- unname(model.matrix(~0+as.factor(trt)))
  des_repr <- unname(model.matrix(~0+as.factor(Rep)))
  one_vecr <- matrix(1, nrow = n, ncol = 1)
  des_fulr <- cbind(one_vecr,des_trr,des_repr)

  ## calculation of c and q matrices
  cmat <- (t(des_trr) %*% des_trr) - ((t(des_trr) %*% des_repr) %*% ginv(t(des_repr) %*% des_repr) %*%
                                        (t(des_repr) %*% des_trr))
  qmat <- (t(des_trr) %*% resp) - ((t(des_trr) %*% des_repr) %*% ginv(t(des_repr) %*% des_repr) %*%
                                     (t(des_repr) %*% resp))
  bmat <- I_mat - des_repr %*% ginv(t(des_repr) %*% des_repr) %*% t(des_repr)

  ## function for generating u matrix
  generate_umatrix <- function(n, t) {
    # Initialize a list to store the matrices
    matrix_list <- list()
    if (t == 1) {
      # Special case when t = 1: simply iterate over each row
      for (i in 1:n) {
        U <- matrix(0, nrow = n, ncol = 1)
        U[i, 1] <- 1  # Place 1 in each row for the single column
        ##storing to a matrix
        matrix_list[[i]] <- U
      }
    } else {
      # Generate all combinations where each column can take any row position (1:n)
      row_combinations <- expand.grid(rep(list(1:n), t))

      # Filter combinations to ensure the first (t-1) columns are in ascending order
      row_combinations <- row_combinations[do.call(order, as.list(row_combinations)), ]

      # Initialize a counter for combinations
      counter <- 1

      # Loop through each valid combination
      for (comb in 1:nrow(row_combinations)) {
        # Initialize an empty matrix U with dimensions n x t
        U <- matrix(0, nrow = n, ncol = t)

        # Place 1's in the matrix for the current combination
        for (j in 1:t) {
          U[row_combinations[comb, j], j] <- 1  # Place 1 in the selected row for each column
        }

        # Store the matrix in the result list
        matrix_list[[counter]] <- U
        counter <- counter + 1
      }
    }
    return(matrix_list)
  }

  ## now generating u matrix for given t
  U_mat <- generate_umatrix(n,noutli)

  ## requirements for cook distance
  S_mat <- bmat %*% des_trr %*% pinv(cmat) %*% t(des_trr) %*% bmat
  V_mat <- bmat - S_mat

  ## for cook distance
  ckd_mat <- matrix(0,nrow = length(U_mat), ncol = 1)

  for (i in 1:length(U_mat)){
    hdel <- ginv(t(U_mat[[i]]) %*% V_mat %*% U_mat[[1]]) %*% t(U_mat[[i]]) %*% V_mat %*% resp
    cknum <- t(hdel) %*% t(U_mat[[i]]) %*% S_mat %*% U_mat[[i]] %*% hdel
    ckden <- (v-1) * var(resp)
    ckdis <- cknum/ckden
    ckd_mat[i,] <-  ckdis
  }

  ## to mark values with threshold values
  ## function to mark
#' Title
#'
#' @param mat
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
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
  cat('------------------------------------------------------------------------
      \nCook s Distance : \n')
  print(ckd_matt)
}
