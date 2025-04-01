#need to create dummy variables in education for logistic regression
#make 1st level of education a reference category
one_hot_encode <- function(df, cat_var) {
  df[[cat_var]] <- droplevels(factor(df[[cat_var]]))  # Ensure it's a factor
  dummies <- model.matrix(~ df[[cat_var]], data = df)[, -1, drop = FALSE]  # Remove first column (reference)
  colnames(dummies) <- sub("df\\[\\[cat_var\\]\\]", cat_var, colnames(dummies))  #  # Find the original column index
  col_index <- which(names(df) == cat_var)
  
  # Create a new dataframe with the dummy columns inserted at the original position
  df_new <- cbind(df[ , 1:(col_index - 1), drop = FALSE],  # Columns before categorical
                  dummies,  # One-hot encoded variables
                  df[ , (col_index + 1):ncol(df), drop = FALSE])  # Columns after
  
  return(df_new)
}


# function to set up Gibbs Sampler for beta using Polya-Gamma augmentation
gibbs_sampler <- function(y, X, n_iter = 5000, b_prior = NULL, B_prior_inv = NULL) {
  
  # Data dimensions
  N <- length(y)  # Number of observations
  P <- ncol(X)    # Number of predictors (including intercept)
  
  # Set priors if not provided
  if (is.null(b_prior)) {
    b_prior <- rep(0, P)  # Mean of Gaussian prior
  }
  if (is.null(B_prior_inv)) {
    B_prior_inv <- diag(1, P)  # Precision matrix 
  }
  
  # Initialize storage for beta samples
  beta_samples <- matrix(0, nrow = n_iter, ncol = P)
  
  # Initial values
  set.seed(123)
  beta_init<-rnorm(P, mean=0, sd=0.01)
  beta <- beta_init  # Start with values close to zero
  
  # Gibbs sampling loop
  for (iter in 1:n_iter) {
    
    # Compute psi = X * beta (log-odds)
    psi <- X %*% beta
    
    # Sample omega from Polya-Gamma distribution
    omega <- rpg(N, rep(1, N), psi)  # PG(1, psi) for logistic regression
    
    # Construct diagonal matrix Omega
    #Omega <- diag(as.vector(omega), N, N)
    
    # Compute kappa
    kappa <- y - 0.5  # Equivalent to (y - n/2) since n = 1
    
    # Compute posterior covariance and mean
    # t(X) %*% (X * omega)
    A       <- t(X) %*%(X*omega)+ B_prior_inv
    R      <- chol(A)
    R.inv   <-backsolve(R, diag(rep(1,ncol(R))))
    V_omega <- R.inv%*%t(R.inv)
    m_omega <- V_omega %*% (t(X) %*% kappa + B_prior_inv %*% b_prior)
    
    # Generate N(0,1) deviates in length of m_omega
    z_i <- rnorm(length(m_omega), 0, 1)
    #apply transformation so beta~MVN(m_omega,V_omega)
    beta <- (R.inv%*%z_i)+m_omega
    
    # Store the sample
    beta_samples[iter, ] <- beta
  }
  
  return(beta_samples)
}

logistic_posterior_predictive_vec <- function(Y, X, y, x_new, theta) {
  # Y: vector of binary outcomes (0,1) for existing observations
  # X: design matrix where each row corresponds to an observation in Y
  # y: reference outcome value (0,1) for computing weights
  # x_new: design vector for new point (for computing weights)
  # theta: MCMC parameter samples (each column is a sample)
  
  n_samples <- ncol(theta)
  n_obs <- length(Y)
  
  # Calculate x_new probabilities once for all samples
  log_probs_new <- drop(1 / (1 + exp(-x_new %*% theta)))
  
  # Pre-allocate matrix for all probabilities
  all_probs <- matrix(0, nrow=n_obs, ncol=n_samples)
  
  if(y == 1) {
    # Normalize weights for y=1
    normalized_weights <- log_probs_new / sum(log_probs_new)
    
    # Compute probabilities for each observation and sample
    for(i in 1:n_obs) {
      all_probs[i,] <- drop(1 / (1 + exp(-X[i,,drop=FALSE] %*% theta)))
    }
    
  } else { # y == 0
    # Normalize weights for y=0
    normalized_weights <- (1 - log_probs_new) / sum(1 - log_probs_new)
    
    # Compute probabilities for each observation and sample
    for(i in 1:n_obs) {
      all_probs[i,] <- drop(1 / (1 + exp(X[i,,drop=FALSE] %*% theta)))
    }
  }
  
  # Apply weights to each column (each sample)
  weighted_probs <- sweep(all_probs, 2, normalized_weights, "*")
  
  # Sum across samples for each observation
  weighted_sums <- rowSums(weighted_probs)
  
  # Final result based on Y values
  results <- ifelse(Y == 1, weighted_sums, 1 - weighted_sums)
  
  return(results)
}

#modify function to output string containing prediction set
full_conformal_classify <- function(Y, X, y_grid, x_new, theta, alpha) {
  # Track which values are accepted
  accepted_0 <- FALSE
  accepted_1 <- FALSE
  
  for (l in 1:length(y_grid)) {
    # Conformity scores on dataset
    sig_1_to_n <- logistic_posterior_predictive_vec(Y = Y, X = X, y = y_grid[[l]], x_new = x_new, theta = theta)
    
    # Conformity score on test point
    sig_n_plus_one <- logistic_posterior_predictive_vec(Y = y_grid[[l]], X = x_new, y = y_grid[[l]], x_new = x_new, theta = theta)
    # Adjusted quantile calculation
    n <- length(Y)
    pi <- (length(which(sig_1_to_n <= sig_n_plus_one)) + 1)/(n+1) 
    # Reject points if pi <= alpha
    if (pi > alpha) {
      # Mark the value as accepted
      if (y_grid[[l]] == 0) {
        accepted_0 <- TRUE
      } else if (y_grid[[l]] == 1) {
        accepted_1 <- TRUE
      }
    }
  }
  
  # Format the output string based on which values were accepted
  if (accepted_0 && accepted_1) {
    return("{0,1}")
  } else if (accepted_0) {
    return("{0}")
  } else if (accepted_1) {
    return("{1}")
  } else {
    return("{}")
  }
}
#parallel functrion computes conformal sets over x_new
conf.sets.parallel <- function(x_new, train.resp, train.data, y.grid, beta, alpha) {
  # Create cluster only once
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export necessary functions
  clusterExport(cl, varlist = c("logistic_posterior_predictive_vec", "full_conformal_classify"), 
                envir = environment())
  
  # Distribute row indices across cores
  row_indices <- 1:nrow(x_new)
  chunks <- split(row_indices, cut(row_indices, num_cores, labels=FALSE))
  
  conf_sets <- foreach(idx = chunks, .combine = c, .packages = c("foreach")) %dopar% {
    # Create storage for results from this chunk
    chunk_results <- numeric(length(idx))
    
    # Process each row in this chunk
    for (i in seq_along(idx)) {
      row_idx <- idx[i]
      x_new_i <- x_new[row_idx, , drop = FALSE]
      
      # Assuming full_conformal_classify returns a single value per row
      chunk_results[i] <- full_conformal_classify(
        Y = train.resp, 
        X = train.data, 
        y_grid = y.grid, 
        x_new = x_new_i,
        theta = beta, 
        alpha = alpha
      )
    }
    
    return(chunk_results)
  }
  
  stopCluster(cl)
  return(conf_sets)
}
#computes Bayesian prediction sets
bayes.pred <- function(x_new, theta, alpha) {
  # Efficiently compute all likelihoods in one operation
  # For each row in x_new, compute probability for each column in theta
  n_rows <- nrow(x_new)
  n_cols <- ncol(theta)
  
  # Pre-allocate matrix to store all probabilities
  all_probs <- matrix(0, nrow=n_rows, ncol=n_cols)
  
  # Compute probabilities for all observations and all theta samples
  for(i in 1:n_rows) {
    linear_pred <- x_new[i,] %*% theta
    all_probs[i,] <- 1 / (1 + exp(-linear_pred))
  }
  
  # Calculate mean probability for each row
  p_means <- rowMeans(all_probs)
  
  # Create result vector
  results <- vector("character", n_rows)
  
  # Assign prediction sets based on conditions
  results[1-p_means >= 1-alpha] <- "{0}"
  results[p_means >= 1-alpha] <- "{1}"
  results[pmax(1-p_means, p_means) <= 1-alpha] <- "{0,1}"
  
  return(results)
}
#optimizes uniformiativ rate for alpha 
plot_uninformative_rates <- function(train.resp, train.data, y.grid, x_new, theta, 
                                     alpha_values = seq(0.01, 0.5, by = 0.05)) {
  
  # Create storage for results
  results <- data.frame(
    alpha = alpha_values,
    empty_count = 0,
    full_count = 0,
    uninformative_rate = 0
  )
  
  # For each alpha value
  for (a_idx in 1:length(alpha_values)) {
    alpha <- alpha_values[a_idx]
    cat(sprintf("Processing alpha = %.2f (%d/%d)\n", alpha, a_idx, length(alpha_values)))
    
    # Use the parallel version to get all confidence sets at once
    store <- conf.sets.parallel(
      x_new = x_new,
      train.resp = train.resp, 
      train.data = train.data, 
      y.grid = y.grid, 
      beta = theta, 
      alpha = alpha
    )
    
    # Count empty and full sets
    empty_count <- sum(store == "{}")
    full_count <- sum(store == "{0,1}")
    uninformative_count <- empty_count + full_count
    
    # Store results
    results$empty_count[a_idx] <- empty_count
    results$full_count[a_idx] <- full_count
    results$uninformative_rate[a_idx] <- uninformative_count/nrow(x_new)
    
    cat(sprintf("  Done: empty=%d, full=%d, uninformative_rate=%.4f\n", 
                empty_count, full_count, uninformative_count/nrow(x_new)))
  }
  
  # Plot the results
  par(mfrow=c(1,2))
  
  # Plot uninformative rate
  plot(results$alpha, results$uninformative_rate, type="o", col="blue", 
       xlab="Alpha", ylab="Uninformative Prediction Rate",
       main="Rate of Uninformative Predictions", ylim=c(0.1,0.9))
  grid()
  
  # Plot breakdown of empty vs full sets
  barplot(t(as.matrix(results[,c("empty_count", "full_count")])), 
          beside=TRUE, names.arg=results$alpha,
          col=c("red", "green"), 
          main="Breakdown of Uninformative Predictions",
          xlab="Alpha", ylab="Count")
  legend("topright", legend=c("Empty sets", "Full sets {0,1}"), 
         fill=c("red", "green"))
  
  # Reset plot settings
  par(mfrow=c(1,1))
  
  return(results)
}
compare_missclassify_rates <- function(train.resp, train.data, y.grid, x_new, y_true, theta, 
                                       alpha_values = seq(0.01, 0.5, by = 0.05)) {
  
  # Create storage for conformal results
  conf_results <- data.frame(
    alpha = alpha_values,
    empty_count = 0,
    full_count = 0,
    single_element_count = 0,
    misclassification_count = 0,
    uninformative_rate = 0,
    misclassification_rate = 0
  )
  
  # Create storage for Bayes results
  bayes_results <- data.frame(
    alpha = alpha_values,
    empty_count = 0,
    full_count = 0,
    single_element_count = 0,
    misclassification_count = 0,
    uninformative_rate = 0,
    misclassification_rate = 0
  )
  
  # For each alpha value
  for (a_idx in 1:length(alpha_values)) {
    alpha <- alpha_values[a_idx]
    cat(sprintf("Processing alpha = %.2f (%d/%d)\n", alpha, a_idx, length(alpha_values)))
    
    # --------- Conformal Method ---------
    cat("  Running conformal method...\n")
    
    # Use the parallel version to get all confidence sets at once
    conf_store <- conf.sets.parallel(
      x_new = x_new,
      train.resp = train.resp, 
      train.data = train.data, 
      y.grid = y.grid, 
      beta = theta, 
      alpha = alpha
    )
    
    # Count empty and full sets for conformal
    conf_empty_count <- sum(conf_store == "{}")
    conf_full_count <- sum(conf_store == "{0,1}")
    conf_uninformative_count <- conf_empty_count + conf_full_count
    
    # Count single-element predictions and misclassifications for conformal
    conf_single_0_indices <- conf_store == "{0}"
    conf_single_1_indices <- conf_store == "{1}"
    conf_single_element_count <- sum(conf_single_0_indices) + sum(conf_single_1_indices)
    
    # Calculate misclassifications for conformal
    conf_misclassifications <- sum(conf_single_0_indices & y_true == 1) + 
      sum(conf_single_1_indices & y_true == 0)
    
    # Store conformal results
    conf_results$empty_count[a_idx] <- conf_empty_count
    conf_results$full_count[a_idx] <- conf_full_count
    conf_results$single_element_count[a_idx] <- conf_single_element_count
    conf_results$misclassification_count[a_idx] <- conf_misclassifications
    conf_results$uninformative_rate[a_idx] <- conf_uninformative_count/nrow(x_new)
    
    # Calculate misclassification rate for conformal
    if (conf_single_element_count > 0) {
      conf_results$misclassification_rate[a_idx] <- conf_misclassifications/conf_single_element_count
    } else {
      conf_results$misclassification_rate[a_idx] <- NA
    }
    
    # --------- Bayesian Method ---------
    cat("  Running Bayesian method...\n")
    
    # Preallocate store for Bayes
    bayes_store <- vector("character", length = nrow(x_new))
    
    # Loop through each row of x_new for Bayes
    for (i in 1:nrow(x_new)) {
      # Pass only the current row to bayes.pred
      bayes_store[i] <- bayes.pred(x_new = x_new[i, , drop=FALSE], theta = theta, alpha = alpha)
      
      # Optional progress indicator 
      if (i %% 50 == 0 || i == nrow(x_new)) {
        cat(sprintf("    Progress: %d/%d rows completed\r", i, nrow(x_new)))
        flush.console()
      }
    }
    
    # Count empty and full sets for Bayes
    bayes_empty_count <- sum(bayes_store == "{}")
    bayes_full_count <- sum(bayes_store == "{0,1}")
    bayes_uninformative_count <- bayes_empty_count + bayes_full_count
    
    # Count single-element predictions and misclassifications for Bayes
    bayes_single_0_indices <- bayes_store == "{0}"
    bayes_single_1_indices <- bayes_store == "{1}"
    bayes_single_element_count <- sum(bayes_single_0_indices) + sum(bayes_single_1_indices)
    
    # Calculate misclassifications for Bayes
    bayes_misclassifications <- sum(bayes_single_0_indices & y_true == 1) + 
      sum(bayes_single_1_indices & y_true == 0)
    
    # Store Bayes results
    bayes_results$empty_count[a_idx] <- bayes_empty_count
    bayes_results$full_count[a_idx] <- bayes_full_count
    bayes_results$single_element_count[a_idx] <- bayes_single_element_count
    bayes_results$misclassification_count[a_idx] <- bayes_misclassifications
    bayes_results$uninformative_rate[a_idx] <- bayes_uninformative_count/nrow(x_new)
    
    # Calculate misclassification rate for Bayes
    if (bayes_single_element_count > 0) {
      bayes_results$misclassification_rate[a_idx] <- bayes_misclassifications/bayes_single_element_count
    } else {
      bayes_results$misclassification_rate[a_idx] <- NA
    }
    
    cat("\n")
  }
  
  # Plot the misclassification rates together
  plot(conf_results$alpha, conf_results$misclassification_rate, type="o", col="blue", 
       xlab="Alpha", ylab="Misclassification Rate",
       main="Comparison of Misclassification Rates",
       ylim=c(0, max(c(conf_results$misclassification_rate, bayes_results$misclassification_rate), na.rm=TRUE) * 1.1))
  lines(bayes_results$alpha, bayes_results$misclassification_rate, type="o", col="red")
  legend("topleft", legend=c("Conformal", "Bayesian"), 
         col=c("blue", "red"), lty=1, pch=1)
  grid()
  
  # Return both results as a list
  return(list(
    conformal = conf_results,
    bayesian = bayes_results
  ))
}
compute_single_element_misclassification <- function(predictions, actuals) {
  # Input validation
  if(length(predictions) != length(actuals)) {
    stop("Predictions and actuals must have the same length")
  }
  
  n <- length(predictions)
  
  # Initialize counters
  single_element_count <- 0
  incorrect_count <- 0
  
  for(i in 1:n) {
    pred <- predictions[i]
    actual <- actuals[i]
    
    # Only consider single element predictions
    if(pred == "{0}" || pred == "{1}") {
      single_element_count <- single_element_count + 1
      
      # Check if prediction is correct
      if((pred == "{0}" && actual != 0) || (pred == "{1}" && actual != 1)) {
        incorrect_count <- incorrect_count + 1
      }
    }
  }
  
  # Calculate misclassification rate for single element predictions
  misclassification_rate <- if(single_element_count > 0) {
    incorrect_count / single_element_count
  } else {
    NA  # No single element predictions to evaluate
  }
  
  return(list(
    single_element_predictions = single_element_count,
    incorrect_predictions = incorrect_count,
    misclassification_rate = misclassification_rate,
    single_element_proportion = single_element_count / n
  ))
}

