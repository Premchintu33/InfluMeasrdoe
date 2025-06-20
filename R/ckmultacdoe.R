#' To identify Outliers in Multi-Response Experiments with Auto-Correlated Errors(1)
#'
#' @param trt Numeric or complex vector containing the treatment levels
#' @param Rep Numeric or complex vector containing the Replication
#' @param resps data frame containing response variables
#' @param noutli a  number indicating the number of outliers
#'
#' @importFrom stats xtabs acf shapiro.test var cov lm pchisq sd IQR model.matrix residuals
#'
#'
#' @return The output contains Multivariate Normality test and BoxM test for
#' Residuals of the of the model and Cook's distance for each treatment or a
#' combination of treatments in Multi - Response Experiments with Auto-Correlated Errors (1)
#'
#' @examples
#' # example code
#' data(ex2)
#' ckmultacdoe(ex3$trt,ex3$rep,ex3[,3:4],2)
#'
#' @export

ckmultacdoe <- function(trt,Rep,resps,noutli){
  argmts <- lapply(as.list(match.call())[-1], eval)
  missing_vals <- sapply(argmts, function(x) any(is.na(x)))
  ## labeling
  resps <- as.matrix(resps)
  n <- length(trt)
  p <- ncol(resps)
  I_mat <- diag(n)
  I_p <- diag(p)
  I_np <- I_p %x% I_mat
  r <- length(unique(Rep)) ## number of replications
  v <- length(unique(trt)) ## number of treatments
  j.vec <-  matrix(1, nrow = n, ncol = 1)
  varcovar <- cov(resps)

  ## checking for any missing values
  if(any(missing_vals)){
    warning("One or more arguments contain missing values.")
  }

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
  p_str <-  generate_contrast_matrix(v)
  P <- (I_p %x% p_str)

  ##design matrices for single response
  des_mtrr <- unname(model.matrix(~0+as.factor(trt)))
  des_repr <- unname(model.matrix(~0+as.factor(Rep)))
  one_vecr <- matrix(1, nrow = n, ncol = 1)
  des_rep1 <- cbind(one_vecr,des_repr)
  des_fulr <- cbind(one_vecr,des_mtrr,des_repr)

  ##Recquired design matrices for multi response
  des_mtrr <- I_p %x% des_mtrr
  des_mrep1 <- I_p %x% des_rep1

  omg_m <- list()

  ## calculation of auto - correlation
  acf_mres <- matrix(0, nrow = length(unique(Rep)), ncol = 1)
  for(j in 1:length(unique(Rep))){
    klm <-  xtabs(resps[,1]~trt+Rep)
    kat=acf(klm[,j], lag.max = NULL, plot = FALSE)
    acf_mres[j] <- kat$acf[2]
  }

  omg_m <- list()
  for (k in 1:nrow(acf_mres)){
    omg_m[[k]] <- matrix(0, nrow = length(unique(trt)), ncol =length(unique(trt)))##no of trt
  }

  for (k in 1:nrow(acf_mres)){
    for (i in 1:length(unique(trt))){
      for (j in 1:length(unique(trt))){
        rho_m <- (1/(1-acf_mres[k]^2))
        if (i==j){
          omg_m[[k]][i,j]=1
        }
        if (i < j){
          omg_m[[k]][i,j] = acf_mres[k]^(j-i)
        }
        if (i >j) {
          omg_m[[k]][i,j] = acf_mres[k]^(i-j)
        }
      }
    }
    omg_m[[k]] <- rho_m * omg_m[[k]]
  }

  ## diagonally combining all the omg matrices
  omg_mfl <- as.matrix(Matrix::bdiag(omg_m))
  OMG <- varcovar %x% omg_mfl


  ## calculation of c and q matrices for  responses
  bmat <- (I_np - des_mrep1 %*% MASS::ginv(t(des_mrep1)  %*% des_mrep1)%*%
              t(des_mrep1))
  cmat <- t(des_mtrr) %*% bmat %*% des_mtrr
  qmat <- t(des_mtrr) %*% bmat %*% as.vector(resps)

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
    trt_combs <- data.frame(0,nrow = length(explored), ncol = 1)
    for (i in seq_along(explored)) {
      itrt_comb <-  data.frame(0,nrow =1 , ncol = delete_count)
      for (j in seq_len(nrow(explored[[i]]))) {
        itrt_comb[,j] <- paste0("(", explored[[i]][j, "row"], ", ", explored[[i]][j, "col"], ") ")
      }
      trt_combs[i,] <- itrt_comb
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

    return(list(trt_combs = trt_combs[,1:delete_count], marked_matrices = marked_matrices))
  }

  ## exp comb output
  datam_fmat <- xtabs(resps[,1]~trt+Rep)

  exp_output <- explore_combinations(matrix(datam_fmat, nrow = v,ncol = r), delete_count = noutli)
  trmac_comb <- exp_output$trt_combs
  mark_acmats <- exp_output$marked_matrices

  ## loop for obtaining u mat
  U_mats <- list()
  for(i in seq_along(exp_output$marked_matrices)){
    u_vec <- as.vector(mark_acmats[[i]])
    U_mats[[i]] <- u_vec
  }

  ##  for cook distance
  ckdm_mat <- matrix(0, nrow = length(U_mats), ncol = 1)

  for (outl in 1:length(U_mats)){
    u_vec <- matrix(1, nrow = p, ncol = 1) %x% U_mats[[outl]]
    W_mat <- bmat - bmat %*% tensr::mhalf(pracma::inv(OMG)) %*% des_mtrr %*% pracma::pinv(cmat) %*% t(des_mtrr) %*% tensr::mhalf(pracma::inv(OMG)) %*% bmat
    uwu <- t(u_vec) %*% tensr::mhalf(pracma::inv(OMG)) %*% W_mat %*% tensr::mhalf(pracma::inv(OMG)) %*% u_vec
    p_thdiff <- (P %*% pracma::pinv(cmat) %*% t(des_mtrr) %*% tensr::mhalf(pracma::inv(OMG)) %*% bmat %*% tensr::mhalf(pracma::inv(OMG)) %*% u_vec %*% MASS::ginv(uwu) %*%
                 t(u_vec) %*% tensr::mhalf(pracma::inv(OMG)) %*% W_mat %*% tensr::mhalf(pracma::inv(OMG)) %*% as.vector(resps))
    pcp <- P %*% cmat %*% t(P)
    numr <- t(p_thdiff) %*% pracma::inv(pcp) %*% p_thdiff
    denom <- (p*(v-1))
    mckd <- numr/denom
    ckdm_mat[outl,] <-  round(mckd,4)
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

  trmac_comb2 <- as.data.frame(trmac_comb)
  for (i in seq_along(trmac_comb2)) {
    colnames(trmac_comb2)[i] <- paste("(trt_rep)", i, sep = "_")
  }
  ckdmac_df <- cbind.data.frame(trmac_comb2,ckdm_matt)

  cat('------------------------------------------------------------------------
      \nCook s Distance : \n')
  print(ckdmac_df)
  cat('\n------------------------------------------------------------------------\n\n')
  cat('------------------------------------------------------------------------
      \nTreatments Suspected as Outliers : \n')
  print(ckdmac_df[ckdmac_df[[noutli+2]]== "*",])
}

ckmultacdoe(ex3$trt,ex3$rep,ex3[,3:4],2)

