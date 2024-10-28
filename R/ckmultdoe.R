#' To Identify Outliers in Multi - Response Experiments
#'
#' @param trt Numeric or complex vector containing the treatment levels
#' @param Rep Numeric or complex vector containing the Replication
#' @param resps data frame containing response variables
#' @param noutli a  number indicating the number of outliers
#'
#' @importFrom MASS ginv
#' @importFrom Matrix bdiag
#' @importFrom pracma inv pinv
#' @importFrom utils combn
#' @importFrom stats xtabs acf shapiro.test var cov lm pchisq sd IQR model.matrix
#' @importFrom tensr mhalf
#'
#' @return The output contains Shapiro-Wilk Normality test and Bartlett test for
#' Residuals of the of the model and Cook's distance for each treatment or a
#' combination of treatments in Multi - Response Experiments
#'
#' @examples
#' data(ex2)
#' attach(ex2)
#' ckmultdoe(ex2$trt,ex2$rep,ex2[,3:4],1)
#'
#' @export

ckmultdoe <- function(trt,Rep,resps,noutli){
  argmts <- lapply(as.list(match.call())[-1], eval)
  missing_vals <- sapply(argmts, function(x) any(is.na(x)))
  ## labeling
  resps <- as.matrix(resps)
  n <- length(trt)
  p <- ncol(resps)
  I_mat <- diag(n)
  I_p <- diag(p)
  v <- length(unique(trt)) ## number of treatments
  j.vec <-  matrix(1, nrow = n, ncol = 1)
  varcovar <- cov(resps)
  Omega <- varcovar %x% I_mat
  ## checking for any missing values
  if(any(missing_vals)){
    warning("One or more arguments contain missing values.")
  }
  ## bartlett code
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
  ## checking for normality and homogeneity of variances
  trt<-as.factor(trt)
  Rep<-factor(Rep)
  n_result <- matrix(0, nrow = ncol(resps), ncol = 1)
  counter <- 1
  for (resp in 1:ncol(resps)){
    model<-lm(as.matrix(resps[,resp])~Rep+trt)
    n_test<-shapiro.test(model$residuals)
    n_result[resp,] <-n_test$p.value
    cat('------------------------------------------------------------------------\n')
    if (n_result[resp,] > 0.05){
      cat('Normality Assumption of response' , resp, ': \n')
      print("At 5% of significance, Normality Assumption is not violated")
    } else {
      cat('Normality Assumption of response' , resp, ': \n')
      print("At 5% of significance, Normality Assumption is  violated")
    }
    cat('\n------------------------------------------------------------------------\n\n')
    cat('------------------------------------------------------------------------\n')
    bartmout <- bartlett(trt,resps[,resp],length(unique(trt)),length(unique(Rep)))
    if (bartmout > 0.05){
      cat('Homogenity of variances of response', resp, ': \n')
      print("At 5% of significance, residuals can  be considered homocedastic!")
    } else {
      cat('Homogenity of variances of response', resp, ': \n')
      print("At 5% of significance, residuals can not be considered homocedastic!")
    }
    cat('\n------------------------------------------------------------------------\n\n')
  }

  ## for P matrix
  fl <-  diag(v) ## for creating p.str matrix
  p_str <- fl[-c(1),]
  P <- (diag(p) %x% p_str)

  ##design matrices
  des_trr <- unname(model.matrix(~0+as.factor(trt)))
  des_repr <- unname(model.matrix(~0+as.factor(Rep)))
  one_vecr <- matrix(1, nrow = n, ncol = 1)
  des_rep1 <- cbind(one_vecr,des_repr)
  des_fulr <- cbind(one_vecr,des_trr,des_repr)

  ## calculation of c and q matrices for multi responses
  bmats <- (I_mat - des_rep1 %*% ginv(t(des_rep1)  %*% des_rep1)%*%
              t(des_rep1))
  cmat <- t(des_trr) %*% bmats %*% des_trr
  qmat <- t(des_trr) %*% bmats %*% resps[,1]

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

  ## generating U matrix for multi responses
  U_mats <- generate_umatrix(n,noutli)

  ## requirements for cook distance
  S_mats <-  bmats %*% des_trr %*% pinv(cmat) %*% t(des_trr) %*% bmats
  V_mats <- bmats - S_mats
  W_lis <- list()

  for (outl in 1:n){
    W_str <- (p_str %*% pinv(cmat) %*% t(des_trr) %*% bmats %*% U_mats[[outl]] %*%
                ginv(t(U_mats[[outl]]) %*% V_mats %*% U_mats[[outl]]) %*% t(U_mats[[outl]]) %*% V_mats)
    W_lis[[outl]] <- W_str
  }

  ckdm_mat <- matrix(0, nrow = (2*n), ncol = 1)
  resps <- as.matrix(resps)
  for(outli in 1:(2*n)){
    index <- ifelse(outli <= n, outli, (outli - 1) %% n + 1)
    ## cook distance
    nume <- t(as.vector(resps)) %*% (inv(varcovar) %x% (t(W_lis[[index]]) %*% (p_str %*% cmat %*% t(p_str)) %*%
                                                          W_lis[[index]])) %*% as.vector(resps)
    denom <- (p*(v-1))

    mckd <- nume/denom
    ckdm_mat[outli,] <-  mckd
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
  fthdm <- IQR(ckdm_mat)
  trshldm <- (7/2) * fthdm
  ckdm_matt <- mark_values_with_asterisk(ckdm_mat,trshldm)
  cat('------------------------------------------------------------------------
      \nCook s Distance : \n')
  print(ckdm_matt)
}
