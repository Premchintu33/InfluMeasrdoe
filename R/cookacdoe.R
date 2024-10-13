ckacdoe <- function(trt,Rep,resp,noutli){
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
  v <- length(unique(trt)) ## number of treatments
  j.vec <-  matrix(1, nrow = n, ncol = 1)

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
    kat=acf(klm[,j], lag = 1, pl=FALSE)
    acf_res[j] <- kat$acf[2]
  }

  omg <- list()
  for (k in 1:nrow(acf_res)){
    omg[[k]] <- matrix(0, nrow = length(unique(trt)), ncol =length(unique(trt)))##no of trt
  }

  for (k in 1:nrow(acf_res)){
    for (i in 1:length(unique(trt))){
      for (j in 1:length(unique(trt))){
        rho <- (1/(1-acf.res[k]^2))
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
  omg_fl <- as.matrix(bdiag(omg))

  ## calculation of c and q matrices
  bmat_ac <- inv(omg_fl) - inv(omg_fl) %*% des_repr %*% ginv(t(des_repr) %*% inv(omg_fl) %*% des_repr) %*% t(des_repr) %*% inv(omg_fl)

  cmat_ac <- t(des_trr) %*% bmat_ac %*% des_trr
  qmat_ac <- t(des_trr) %*% bmat_ac %*% resp

  ## requirements for cook distance
  Vmat_ac <- omg_fl %*% (bmat_ac - bmat_ac %*% des_trr %*% ginv(cmat_ac) %*% t(des_trr) %*% bmat_ac)
  Hstr <- inv(omg_fl) %*% Vmat_ac

  ## for cook distance
  ckd_acmat <- matrix(0,nrow =n, ncol = 1)

  # Function to generate all possible outlier vectors with specified number of outliers
  generate_outvec <- function(n, num_outliers) {
    vect_list <- list()
    # Get all possible combinations of positions for the given number of outliers
    outlier_combinations <- combn(1:n, num_outliers)

    # Loop through each combination of outlier positions
    for (i in 1:ncol(outlier_combinations)) {
      # Initialize a vector of zeros
      U <- rep(0, n)

      # Get the current combination of outlier positions
      outliers <- outlier_combinations[, i]

      # Set the corresponding positions in U to 1
      U[outliers] <- 1
      vect_list[[i]] <- U
    }
    return(vect_list)
  }

  u_lis <- generate_outvec(n,noutli)

  ## loop for cook's distance
  for(outl in 1:n){
    omg_str_nu <- t(as.vector(u_lis[[outl]])) %*% inv(omg_fl) %*% Vmat_ac %*% resp
    omg_str_de <- sd(resp) * sqrt(t(as.vector(u_lis[[outl]])) %*% Vmat_ac %*% as.vector(u_lis[[outl]]))
    omg_str <- omg_str_nu/omg_str_de
    ckd_ac_fr <- (bmat_ac[outl,outl] - Hstr[outl,outl])/Hstr[outl,outl]
    ckd_ac_sc <- (omg_str^2)/(v-1)
    ckd_ac <- ckd_ac_fr * ckd_ac_sc
    ckd_acmat[outl,] <- ckd_ac
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
  cat('------------------------------------------------------------------------
      \nCook s Distance : \n')
  print(ckd_acmatt)
}
