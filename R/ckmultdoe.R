#' To Identify Outliers in Multi - Response Experiments
#'
#' @param trt Numeric or complex vector containing the treatment levels
#' @param Rep Numeric or complex vector containing the Replication
#' @param resps data frame containing response variables
#' @param noutli a  number indicating the number of outliers
#'
#' @importFrom stats xtabs acf shapiro.test var cov lm pchisq sd IQR model.matrix residuals
#'
#' @return The output contains Multivarite Normality test and BoxM test for
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
  r <- length(unique(Rep)) ## number of replications
  v <- length(unique(trt)) ## number of treatments
  j.vec <-  matrix(1, nrow = n, ncol = 1)
  varcovar <- cov(resps)
  Omega <- varcovar %x% I_mat
  ## checking for any missing values
  if(any(missing_vals)){
    warning("One or more arguments contain missing values.")
  }

  ## checking for normality and homogeneity of variances
  trt<-as.factor(trt)
  Rep<-factor(Rep)

  ## model to extract residuals
  trt<-as.factor(trt)
  Rep<-factor(Rep)

  model_fres <- lm(resps ~ trt+Rep)
  resi <- residuals(model_fres) ## for extracting residuals

  mvnres <- MVN::mvn(resi)## multivariate normality test

  cat('------------------------------------------------------------------------\n')
  if (mvnres$multivariateNormality[,3]> 0.05){
    cat('Normality Assumption of response: \n')
    print("At 5% of significance, Normality Assumption is not violated")
    print(mvnres$multivariateNormality)
  } else {
    cat('Normality Assumption \n'):
      print("At 5% of significance, Normality Assumption is  violated")
    print(mvnres$multivariateNormality)
  }
  cat('\n------------------------------------------------------------------------\n\n')
  cat('------------------------------------------------------------------------\n')

  boxres <- heplots::boxM(resi, trt) ## multivariate homogenity of var covar test

  if (boxres$p.value > 0.05){
    cat('Homogenity of variances of response: \n')
    print("At 5% of significance, residuals can  be considered homocedastic!")
    print(boxres)
  } else {
    cat('Homogenity of variances of response: \n')
    print("At 5% of significance, residuals can not be considered homocedastic!")
    print(boxres)
  }
  cat('\n------------------------------------------------------------------------\n\n')

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
  p_str <- generate_contrast_matrix(v)
  P <- (diag(p) %x% p_str)

  ##design matrices
  des_trr <- unname(model.matrix(~0+as.factor(trt)))
  des_repr <- unname(model.matrix(~0+as.factor(Rep)))
  one_vecr <- matrix(1, nrow = n, ncol = 1)
  des_rep1 <- cbind(one_vecr,des_repr)
  des_fulr <- cbind(one_vecr,des_trr,des_repr)

  ## calculation of c and q matrices for single responses
  bmats <- (I_mat - des_rep1 %*% MASS::ginv(t(des_rep1)  %*% des_rep1)%*%
              t(des_rep1))
  cmat <- t(des_trr) %*% bmats %*% des_trr
  qmat <- t(des_trr) %*% bmats %*% resps[,1]

  ## function for trt combinations and u_matrix
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
    trt_combs <-  data.frame(matrix(nrow = length(explored), ncol = delete_count))
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

  ## exp combinations output
  datam_fmat <- xtabs(resps[,1]~trt+Rep)

  exp_moutput <- explore_combinations(matrix(datam_fmat, nrow = v,ncol = r), delete_count = noutli)
  tr_mcomb <- exp_moutput$trt_combs
  mark_mmats <- exp_moutput$marked_matrices

  ## generating U matrix for multi responses
  U_mats <- list()
  for(i in seq_along(mark_mmats)){
    u_ini_mmat <- matrix(0,nrow = n, ncol = noutli)
    for(j in 1:noutli){
      u_vec <- as.vector(mark_mmats[[i]][[j]])
      u_ini_mmat[,j] <- u_vec
    }
    U_mats[[i]] <- u_ini_mmat
  }

  ## requirements for cook distance
  S_mats <-  bmats %*% des_trr %*% pracma::pinv(cmat) %*% t(des_trr) %*% bmats
  V_mats <- bmats - S_mats
  W_lis <- list()

  for (outl in 1:length(U_mats)){
    W_str <- (p_str %*% pracma::pinv(cmat) %*% t(des_trr) %*% bmats %*% U_mats[[outl]] %*%
                MASS::ginv(t(U_mats[[outl]]) %*% V_mats %*% U_mats[[outl]]) %*% t(U_mats[[outl]]) %*% V_mats)
    W_lis[[outl]] <- W_str
  }

  ckdm_mat <- matrix(0, nrow = length(U_mats), ncol = 1)
  resps <- as.matrix(resps)
  for(outli in 1:length(U_mats)){
    index <- ifelse(outli <= n, outli, (outli - 1) %% n + 1)
    ## cook distance
    nume <- t(as.vector(resps)) %*% (pracma::inv(varcovar) %x% (t(W_lis[[index]]) %*% (p_str %*% cmat %*% t(p_str)) %*%
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

  tr_mcomb2 <- as.data.frame(tr_mcomb)
  for (i in seq_along(tr_mcomb2)) {
    colnames(tr_mcomb2)[i] <- paste("(trt_rep)", i, sep = "_")
  }
  ckdm_df2 <- cbind.data.frame(tr_mcomb2,ckdm_matt)

  cat('------------------------------------------------------------------------
      \nCook s Distance : \n')
  print(ckdm_df2)
}
