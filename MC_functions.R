########################## Perturbation ########################################


perturb_x <- function(df, Xvar, x, delta = 0.01)
  # Determines the empirical CDF of a unit's value (x) given the variable (Xvar) 
  # and kerncel, select a symmetrical interval around the proportion of units with a value equal
  # or less than the current value, mapped to a quantile range, which is then used
  # to build the range of a uniform distribution, from which a value is drawn.
{
  x <- as.numeric(x)   
  
  # Specifically for A units with delta=0 for which no B-auxiliary values are merged 
  if (is.na(x)) return(x)
  
  # Proportion of values smaller or equal to current value 
  prop_smaller <- mean(df[[Xvar]] <= x, na.rm = TRUE)
  
  # Situation where x is close to 1 
  if (prop_smaller == 1 || prop_smaller > 1 - delta) ranges <- c(prop_smaller - delta, prop_smaller)
  # Situation where x is close to 0
  else if (prop_smaller == 0 || prop_smaller < delta) ranges <- c(prop_smaller, prop_smaller + delta)
  # Add perturbation to create a range
  else ranges <- c(prop_smaller - delta, prop_smaller + delta)
  
  # Map range to obtain values from quantile function
  limits <- quantile(df[[Xvar]], probs = ranges, na.rm = TRUE)
  
  # Mutate x to new x obtained from uniform distribution with ranges equal
  # to quantile limits
  new_x <- round(runif(1, min = limits[1], max = limits[2]), 0)
  
  return(new_x)
}





########################## PCA #################################################

PCA_B <- function(B, auxiliary)
  # Performs PCA on auxiliaries of B, with parallel analysis to determine the
  # optimal number of components. The eigenvalues are later used to create components
  # from auxiliaries in A.
{
  scaling <- function(df)
    # Standardizes auxiliary variables and returns a matrix.
  {
    mu <- colMeans(df[,auxiliary])   # mean
    Sigma <- apply(df[,auxiliary], 2, sd)  # stdev
    Xscaled <- scale(df[,auxiliary], center = mu, scale = Sigma) # standardize
    
    Xscaled %<>% replace(is.nan(.), 0) # replace NaN's
    
    # if a column has 0 st.dev, make sure to exclude it (also in A)
    zero_sd_idx <- which(apply(Xscaled, 2, function(x) sd(x, na.rm=T) == 0))
    if (length(zero_sd_idx) > 0) Xscaled <- Xscaled[, -zero_sd_idx, drop=F]
    
    return(list(X   = Xscaled,
                idx = zero_sd_idx))
  }
  
  # Get standardized B and potential removals
  B_scaled  <- scaling(B)$X
  zero_idx  <- scaling(B)$idx
  
  # Eigendecomposition on correlation matrix of B
  eig <- eigen(cor(B_scaled))
  
  # Parallel analysis
  ap <- parallel(subject = nrow(B), var = ncol(B_scaled))
  nS <- nScree(x = eig$values, aparallel = ap$eigen$qevpea)
  all_comp <- as.matrix(table(unlist(nS$Components)))
  # If across all methods there is no consensus on the optimal number of components,
  # just select the smallest one, else just select based on majority vote.
  if (all(all_comp == all_comp[1])) { 
    num_components <- as.numeric(rownames(all_comp))[which.min(as.numeric(rownames(all_comp)))]
  } else {
    num_components <- as.numeric(rownames(all_comp)[which(all_comp == max(all_comp))])
  }
  
  # Eigenvalues
  L <- as.matrix(eig$vectors[ , 1:num_components, drop = FALSE])
  
  components <- function(df_scaled, df)
    # Replace auxiliaries in B with principal components.
  {
    PC_X <- df_scaled %*% L 
    colnames(PC_X) <- paste0("PC", 1:ncol(PC_X))
    
    pca_df <- df %>% select(-all_of(auxiliary)) 
    pca_df <- cbind(pca_df, PC_X)
    
    return(pca_df)
  }
  # Sample B with principal components
  B_pca  <- components(B_scaled, B)
  
  return(list(B  = B_pca,
              L  = L,
              zero_idx = zero_idx))
}

PCA_A <- function(A, loadings, zero_idx, auxiliary)
  # Uses eigenvalues from PCA on B to transform auxiliaries in A to principal
  # components. 
{
  scaling <- function(df, zero_sd_idx=zero_idx)
    # Standardizes auxiliary variables and returns a matrix.
  {
    mu <- colMeans(df[,auxiliary])   # mean
    Sigma <- apply(df[,auxiliary], 2, sd)  # stdev
    Xscaled <- scale(df[,auxiliary], center = mu, scale = Sigma) # standardize
    
    Xscaled %<>% replace(is.nan(.), 0) # replace NaN's
    
    # If there are any auxiliaries removed from B, also remove them from A
    if (length(zero_sd_idx) > 0) Xscaled <- Xscaled[, -zero_sd_idx, drop=F]
    
    return(Xscaled)
  }
  
  # Get standardized A
  A_scaled  <- scaling(A)
  
  components <- function(df_scaled, df, L=loadings)
    # Replace auxiliaries in A with principal components.
  {
    PC_X <- df_scaled %*% L 
    colnames(PC_X) <- paste0("PC", 1:ncol(PC_X))
    
    pca_df <- df %>% select(-all_of(auxiliary)) 
    pca_df <- cbind(pca_df, PC_X)
    
    return(pca_df)
  }
  
  # Sample A with principal components
  A_pca  <- components(A_scaled, A)
  
  return(A_pca)
}

#################### Pseudo-weighting (Liu et al. 2023) ########################

pseudo_weights <- function(A, B, auxiliary, PCA = FALSE)
  # Computes pseudo-weights for units in B to correct for selection bias by 
  # propensity score weighting using the auxiliaries from A and B, excluding overlap.
  # Then, population totals of auxiliaries can be estimated by applying the 
  # Horvitz-Thompson estimator using the pseudo-weights.
{
  if (PCA) auxiliary <- grep("PC", colnames(A), value=T)   # Determine names of auxiliaries
  
  overlap_ID <- A$BE_ID[A$delta == 1]     # ID of overlapping units
  
  combined_df <- rbind.fill(A, B)   # Combine A and B fully
  
  non_overlap_df <- combined_df %>% filter(!BE_ID %in% overlap_ID)   # Exclude the overlap
  
  # Define logistic regression model (Z=1 if a unit is in B, Z=0 if in A)
  logreg_model <- as.formula(paste("Z ~", paste(auxiliary, collapse=" + ")))
  # Fit logistic regression model on non-overlap data
  model0 <- glm(logreg_model, family=binomial, data=non_overlap_df)
  
  B$O <- exp(predict(model0, newdata=B, type="link"))   # Estimate odds that a unit is in B based on auxiliaries  
  
  B$w_i <- 1 + (B$d_i - 1) / B$O   # Compute pseudo-weights
  
  B %<>% mutate(across(w_i, ~ ifelse(is.nan(.), 1, .))) # Transform NaN's and Inf's to just 1 (simple solution)
  
  B$w_i <- pmin(pmax(B$w_i, 0.01), 100)   # Determine a reasonable range for the pseudo-weights (can explode)
  
  N_TB <- colSums(sweep(as.matrix(B[,auxiliary]), 1, B$w_i, "*"))  # Obtain estimation of auxiliary totals
  
  return(list(N_TB = setNames(N_TB, auxiliary),
              B = B))
}


################### Calibration weighting (Kim & Tam 2021) #####################


calibration <- function(A, B, U, N_TB, KT_error = F, PCA = F)
  # A sophisticated function that handles two diverse approaches at the same time:
  # 1) KT method with and without a measurement error correction model.
  # 2) GREG estimator (Deville & Särndal 1992) using the estimated population
  # totals of the auxiliaries from the pseudo-weighting function. 
{ 
  # Numerically order ID variable
  B <- B[order(B$BE_ID), ]
  A <- A[order(A$BE_ID), ]
  
  A_KT  <- A  # A used for KT method
  A_Liu <- A  # A used for GREG method
  
  X_vars <- names(N_TB)  # Get names of auxiliaries
  
  #----------------------------#
  ## KT method ##
  #----------------------------#
  
  # Delta=1: calibration weighting is similar to GREG methodology
  
  delta1_idx <- A_KT$delta == 1   # Get delta=1 indices in A
  
  # Determine which units overlap
  match_idx <- match(A_KT$BE_ID[delta1_idx], B$BE_ID) 
  delta1_B <- B[match_idx, ]  
  stopifnot(all(delta1_B$BE_ID == A_KT$BE_ID[delta1_idx])) # Raise error if there is discrepancy
  
  if (KT_error) 
  { # The measurement error model involves using the auxiliary values from B instead
    # of A to compute (X^TX)^-1 least squares Gram matrix.
    A_delta1 <- as.matrix(cbind(A_KT$delta[delta1_idx, drop=F], 
                                delta1_B[, X_vars, drop=F]))
  }
  else 
  { # Without correction for errors, the auxiliary values from A are used.
    A_delta1 <- as.matrix(cbind(A_KT$delta, A_KT[,X_vars]))[delta1_idx, , drop=F]
  }
  
  class(A_delta1) <- "numeric"
  # The sigma Gram matrix is computed as: d_i * X * X^T
  weights_delta1 <- A_KT$d_i[delta1_idx]
  X_weighted <- sweep(A_delta1, 1, weights_delta1, "*")  
  sigma_mat <- crossprod(X_weighted, A_delta1)
  sigma_hat <- tryCatch({ # The Gram matrix can be (close to) singular, in which case
    # it is better to opt for the Moore-Penrose (generalized) inverse.
    solve(sigma_mat)  
  }, error = function(e) {
    ginv(sigma_mat)
  })
  
  # The weights for delta=1 are calibrated to B-auxiliary totals
  B_mat <- as.matrix(cbind(B$delta, B[,X_vars]))
  class(B_mat) <- "numeric"
  Xt <- t(as.matrix(colSums(B_mat)))
  
  # Obtain the calibration weights
  A_KT$w_i[delta1_idx] <- weights_delta1 * Xt %*% sigma_hat %*% t(A_delta1)
  
  
  
  # Delta=0: weights are obtained by d_i * (Nc_hat / Nc)
  
  delta0_idx <- A_KT$delta == 0  # Get delta=0 indices in A
  weights_delta0 <- A_KT$d_i[delta0_idx]  # Get sampling weights
  Nc_hat <- sum(weights_delta0)  # Get estimate of missing data size (Nc_hat)
  # Calculate Nc = N - Nb
  N      <- nrow(U)
  Nb     <- nrow(B)
  Nc     <- N - Nb
  # Obtain the calibration weights
  A_KT$w_i[delta0_idx] <- weights_delta0 * (Nc / Nc_hat)
  
  
  
  #----------------------------#
  ## Liu + GREG method ##
  #----------------------------#
  
  # Obtain population auxiliary totals
  N_TB <- N_TB * (N / sum(B$w_i))
  Tx    <- matrix(c(N, N_TB), ncol=1)
  
  # Compute the sigma Gram matrix
  Xi_t <- as.matrix(cbind(1, A_Liu[,X_vars]))
  Xi <- t(Xi_t)
  di_Xi <- sweep(Xi, 2, matrix(A_Liu$d_i, nrow=1), "*")
  sigma_mat <- di_Xi %*% Xi_t
  sigma_hat <- tryCatch({
    solve(sigma_mat)
  }, error = function(e) {
    ginv(sigma_mat)
  })
  
  # Multiply the auxiliary vectors by their design weights
  weighted_X <- sweep(Xi_t, 1, A_Liu$d_i, "*")

  # Compute calibration weights
  Wi <- weighted_X %*% sigma_hat %*% Tx
  A_Liu$w_i <- as.vector(Wi)

  
  return(list(KT = A_KT,
              Liu = A_Liu))
}


##################### Population total estimates ###############################


regression <- function(A_Liu, A_KT, A_KT_corr, var_selection, PCA = FALSE)
{
  
  # Liu + GREG method
  Treg_A_Liu_Wi <- matrix(A_Liu$w_i, nrow=1) %*% 
    as.matrix(A_Liu[ ,var_selection])
  
  # KT method (without measurement error correction)
  Treg_A_KT_Wi <- matrix(A_KT$w_i, nrow=1) %*% 
    as.matrix(A_KT[,var_selection])
  
  # KT method (with measurement error correction)
  Treg_A_KT_Wi_corr <- matrix(A_KT_corr$w_i, nrow=1) %*% 
    as.matrix(A_KT_corr[,var_selection])
  
  if (!PCA) 
  {
    # HT estimator
    Treg_A_HT <- matrix(A_Liu$d_i, nrow=1) %*% 
      as.matrix(A_Liu[,var_selection])
    
    return(list(Liu     = Treg_A_Liu_Wi,
                KT      = Treg_A_KT_Wi,
                KT_corr = Treg_A_KT_Wi_corr,
                HT      = Treg_A_HT))
  }
  
  else
  {
    return(list(Liu_PCA     = Treg_A_Liu_Wi,
                KT_PCA      = Treg_A_KT_Wi,
                KT_corr_PCA = Treg_A_KT_Wi_corr))
  }
}


################ Relative Bias & Coefficient of Variation #####################


get_mean_mat <- function(method, N)
  # Computes average across all simulations. 
  ## E[ \hat{Y} ]
{
  estimates_sum <- purrr::transpose(method) %>% map(reduce, `+`) # Sum up all estimates of a method              
  sum_mat <- do.call(rbind, estimates_sum) # Merge them in a matrix
  rownames(sum_mat) <- names(estimates_mat) # Fix names
  mean_mat <- sum_mat / N   # Compute average over all iterations (N)
  
  return(mean_mat)
}

relative_bias <- function(T_y, method, N)
  # Computes bias relative to true population totals ( ~ accuracy).
  ## [ E[ \hat{Y} ] - Y ] / Y
{
  mean_mat <- get_mean_mat(method, N)  # Get mean estimates across simulation
  relbias <- (mean_mat - T_y) / T_y  # Compute relative bias
  rownames(relbias) <- sapply(method, names)[,1]  # Use the correct method name
  
  return(relbias)
}

coef_of_variation <- function(T_y, method, N)
  # Computes the coefficient of variation of estimates ( ~ precision).
  ## sqrt[ E[ (\hat{Y} - E[\hat{Y}])^2 ] ] / Y
{
  reverse_list <- purrr::transpose(method)   # Flip the list of estimates 
  mean_mat <- get_mean_mat(method, N)  # Get mean estimates across simulation 
  
  variance <- map(1:N, function(j) {  # Map over each simulation/estimate
    df_list <- map(reverse_list, function(x) x[[j]]) # Take simulation j with its estimates
    df <- as.matrix(do.call(rbind, df_list))  # Make it a df (or matrix)
    mat <- (df - mean_mat)^2  # Compute squared difference with mean
    rownames(mat) <- names(reverse_list)  # Fix names
    mat  # Return the squared difference
  }) |> 
    reduce(`+`) |>  # Sum all the squared differences
    divide_by(N - 1)  # Divide by N-1 to obtain the variance f estimates
  
  # Take the square root and divide by true totals to obtain CV
  coef_of_var <-  sqrt(variance) / T_y  
  
  return(coef_of_var)
}




############################### MONTE CARLO ####################################


monte_carlo <- function(population_df, B, N, T_y_mat, T_y_list, auxiliary, var_selection, 
                        A_prop, B_prop, A_slicing = F, B_slicing = F)
  # Monte Carlo combining all individual functions, simulating sample A N-times. 
{
  results <- list()   # list to save results
  
  if (B_slicing) 
  {   # Scenario 5: shrink sample B
    B %<>%   # Randomly remove units per kerncel to achieve the desired B_prop.
      group_by(KERNCEL) %>% 
      slice_sample(prop = B_prop) %>% 
      ungroup()
  }
  
  # Scenario 4: we can scale down A by A_prop by scaling down the sampling
  # weights immediately within the population 
  if (A_slicing) population_df %<>% mutate(d_i = d_i * (1 / A_prop)) 
  
  U_list  <- split(population_df, population_df$KERNCEL)  # Population split by kerncel
  B_list  <- split(B, B$KERNCEL)  # B split by kerncel
  
  
  #----------------------------#
  ## PCA on B-auxiliaries ##  
  #----------------------------#
  pca_B <- map(B_list, ~ PCA_B(.x, auxiliary))  
  B_pca <- lapply(pca_B, function(x) x$B)
  L_pca <- lapply(pca_B, function(x) x$L)
  zero_pca <- lapply(pca_B, function(x) x$zero_idx)
  #----------------------------#
  
  
  for (iter in 1:N)
  {
    if (iter %% 10 == 0) cat("Iteration", paste0(iter, "/", N, "\n")) # Keep track of iterations
    
    #----------------------------#
    ## Sample A ##
    #----------------------------#
    # Take units randomly from population with inclusion probability = 1 / sampling weight,
    # stratified by kerncel.
    group_idx <- population_df %>% 
      group_by(KERNCEL) %>% 
      group_map(~ {
        pik <- 1 / .x$d_i
        .x$BE_ID[UPrandomsystematic(pik) == 1]
      })
    
    sel_units <- do.call(c, group_idx)  
    A <- population_df %>% filter(BE_ID %in% sel_units) 
    
    if (B_slicing) 
    {  # In case of Scenario 5, we need to redefine indicator variable delta
      A %<>% mutate(delta = ifelse(BE_ID %in% B$BE_ID, 1, 0))
    }
    
    A_list  <- split(A, A$KERNCEL) # B split by kerncel
    #----------------------------#



    #----------------------------#
    ## PCA on A-auxiliaries ##
    #----------------------------#
    A_pca <- pmap(list(A_list, L_pca, zero_pca), ~ PCA_A(..1, ..2, ..3, auxiliary))
    #----------------------------#


    #----------------------------#
    ## Pseudo-weighting ##
    #----------------------------#
    # no PCA
    PW              <- map2(A_list, B_list, ~ pseudo_weights(.x, .y, auxiliary))
    N_TB            <- lapply(PW, function(x) x$N_TB)
    B_pseudo_list   <- lapply(PW, function(x) x$B)

    # PCA
    PW_pca              <- map2(A_pca, B_pca, ~ pseudo_weights(.x, .y, auxiliary, PCA=T))
    N_TB_pca            <- lapply(PW_pca, function(x) x$N_TB)
    B_pseudo_pca_list   <- lapply(PW_pca, function(x) x$B)




    #----------------------------#
    ## Calibration weighting
    #----------------------------#
    # no PCA, without measurement error correction + LiuGREG
    cal_clean <- pmap(list(A_list, B_pseudo_list, U_list, N_TB),
                      ~ calibration(..1, ..2, ..3, ..4))

    # no PCA, with measurement error correction + LiuGREG
    cal_error <- pmap(list(A_list, B_pseudo_list, U_list, N_TB),
                      ~ calibration(..1, ..2, ..3, ..4, KT_error=T))

    A_KT      <- lapply(cal_clean, function(x) x$KT)
    A_KT_corr <- lapply(cal_error, function(x) x$KT)
    A_Liu     <- lapply(cal_clean, function(x) x$Liu)
    # A_Liu     <- lapply(cal_error, function(x) x$Liu) similar


    # PCA, without measurement error correction + LiuGREG
    cal_clean_pca <- pmap(list(A_pca, B_pseudo_pca_list, U_list, N_TB_pca),
                          ~ calibration(..1, ..2, ..3, ..4, PCA=T))

    # PCA, with measurement error correction + LiuGREG
    cal_error_pca <- pmap(list(A_pca, B_pseudo_pca_list, U_list, N_TB_pca),
                          ~ calibration(..1, ..2, ..3, ..4, KT_error=T, PCA=T))

    A_KT_pca      <- lapply(cal_clean_pca, function(x) x$KT)
    A_KT_corr_pca <- lapply(cal_error_pca, function(x) x$KT)
    A_Liu_pca     <- lapply(cal_clean_pca, function(x) x$Liu)
    # A_Liu_pca     <- lapply(cal_error_pca, function(x) x$Liu) similar



    #----------------------------#
    ## Population total estimates ##
    #----------------------------#
    # No PCA
    estimates     <- pmap(list(A_Liu, A_KT, A_KT_corr),
                          ~ regression(..1, ..2, ..3, var_selection))

    # PCA
    estimates_pca <- pmap(list(A_Liu_pca, A_KT_pca, A_KT_corr_pca),
                          ~ regression(..1, ..2, ..3, var_selection, PCA = T))

    estimates_full <- Map(c, estimates, estimates_pca)  # Merge them into one list

    results[[iter]] <- Map(c, estimates_full, T_y_list)  # Also merge true population totals
  }


  method_names <- sapply(estimates_full, names)[,1] # Obtain names of all methods
  N <- length(results)  # In case of iterations being skipped


  #----------------------------#
  ## Relative bias ##
  #----------------------------#
  # 1) Map over all methods (map(method_names))
  # 2) Map over each kerncel (map(results))
  # 3) Map over each Y variable within kerncel (map(results ~ map(y[[j]])))
  rel_bias <- map(method_names, function(j) {
    method_results <- map(results, function(x) map(x, function(y) y[[j]]))
    relative_bias(T_y_mat, method_results, N)
  })



  #----------------------------#
  ## Coefficient of variation ##
  #----------------------------#
  # 1) Map over all methods (map(method_names))
  # 2) Map over each kerncel (map(results))
  # 3) Map over each Y variable within kerncel (map(results ~ map(y[[j]])))
  coef_var <- map(method_names, function(j) {
    method_results <- map(results, function(x) map(x, function(y) y[[j]]))
    coef_of_variation(T_y_mat, method_results, N)
  })



  names(rel_bias) <- names(coef_var) <- method_names # Assign correct method names


  #----------------------------#
  ## Relative Mean Squared Error (R-MSE) ##
  #----------------------------#
  MSE <- mapply("+",                                # Sum up:
                lapply(rel_bias, function(x) x^2),  # Squared relative bias
                lapply(coef_var, function(x) x^2),  # Squared coef. of var.
                SIMPLIFY = F)




  #----------------------------#
  ## Relative Root Mean Squared Error (RRMSE) ##
  #----------------------------#
  RRMSE <- lapply(MSE, sqrt)  # Take the root of R-MSE


  return(list(results     = results,
              rel_bias    = rel_bias,
              coef_var    = coef_var,
              RRMSE       = RRMSE))
}