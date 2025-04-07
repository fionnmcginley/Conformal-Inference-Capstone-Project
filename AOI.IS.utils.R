#Gibbs sampler
bayes.SLR.gibbs <- function(n = 100, x, y, alpha_prior_mean = 0, alpha_prior_sd = 1, 
                                     beta_prior_mean = 0, beta_prior_sd = 1, 
                                     iter = 8000, seed = 4, burnin = 4000, sigma) {
  
  set.seed(seed)  # Set seed for reproducibility
  
  # Initialize storage for samples - only need space for post-burnin samples
  alpha_samples <- numeric(iter)
  beta_samples <- numeric(iter)
  
  # Initial values for alpha and beta 
  # Starting closer to true values can reduce burnin needed
  a0 <- rnorm(1, mean = 0, sd = 1)
  b0 <- rnorm(1, mean = 0, sd = 1)
  
  # Gibbs sampling - use separate counter for storing samples
  sample_idx <- 0
  
  # Run the sampler for burnin + iter iterations
  for (i in 1:(burnin + iter)) {
    # Update alpha
    alpha_posterior_var <- 1 / (1 / alpha_prior_sd^2 + n / sigma^2)
    alpha_posterior_mean <- alpha_posterior_var * (
      alpha_prior_mean / alpha_prior_sd^2 + sum(y - b0 * x) / sigma^2
    )
    a0 <- rnorm(1, mean = alpha_posterior_mean, sd = sqrt(alpha_posterior_var))
    
    # Update beta
    beta_posterior_var <- 1 / (1 / beta_prior_sd^2 + sum(x^2) / sigma^2)
    beta_posterior_mean <- beta_posterior_var * (
      beta_prior_mean / beta_prior_sd^2 + sum((y - a0) * x) / sigma^2
    )
    b0 <- rnorm(1, mean = beta_posterior_mean, sd = sqrt(beta_posterior_var))
    
    # Store samples only after burnin
    if (i > burnin) {
      sample_idx <- sample_idx + 1
      alpha_samples[sample_idx] <- a0
      beta_samples[sample_idx] <- b0
    }
  }
  
  # Summarize posterior samples
  alpha_mean <- mean(alpha_samples)
  beta_mean <- mean(beta_samples)
  
  # Add true values to the output for comparison
  return(list(alpha_samples = alpha_samples, 
              beta_samples = beta_samples, 
              alpha_mean = alpha_mean, 
              beta_mean = beta_mean))
}


#function to calculate parameters for analytical posterior distribution
# inputs are prior mean and prior covariance matrix 
posterior.params <- function(X,Y,pri.mu,pri.sigma){
  z<-(1/sigma^2)*t(X)%*%X+solve(pri.sigma)
  theta_mu<-solve(z)%*%(solve(pri.sigma)%*%pri.mu + (1/sigma^2)*t(X)%*%Y)
  list(theta_mu=theta_mu, theta_var_mat=solve(z))
}
#function to plot Gibbs samples
plot.samples <- function(al.samples,be.samples, case="Null", post.params){
  
  # Plot the posterior distributions
  hist(al.samples, breaks=100, main =paste("Posterior of alpha (", case, ")"), xlab = "alpha", freq = F, xlim=c(min(al.samples), max(max(al.samples),2)))
# Create a sequence of x values for the normal distribution curve
a_x_values <- seq(min(al.samples), max(max(al.samples),2), length = 100)

# Calculate the marginal density for those x values
alpha_normal_density <- dnorm(a_x_values, mean = post.params$theta_mu[1], sd = sqrt(post.params$theta_var_mat[1,1]))

# Add the normal distribution density curve to the histogram
lines(a_x_values, alpha_normal_density, col = cbPalette[1], lwd = 2)
# Add a vertical line for the true alpha value
abline(v = 2, col = cbPalette[8], lwd = 2, lty = 2)  # True alpha value

# Plot the  2nd posterior distribution
hist(be.samples, breaks=100, main = paste("Posterior of beta (", case, ")"), xlab = "beta", freq = F,xlim=c(min(be.samples), max(max(be.samples),3)))
# Create a sequence of x values for the normal distribution curve
b_x_values <- seq(min(be.samples), max(max(al.samples),3), length = 100)

# Calculate marginal density for those x values
beta_normal_density <- dnorm(b_x_values, mean = post.params$theta_mu[2], sd = sqrt(post.params$theta_var_mat[2,2]))

# Add the normal distribution density curve to the histogram
lines(b_x_values, beta_normal_density, col = cbPalette[2], lwd = 2)
abline(v = 3, col = cbPalette[8], lwd = 2, lty = 2)  # True beta value
}


#Function to compute upper and lower bounds of conformal interval
#uses parallel computation to loop over new x values
conformal_bounds_single_parallel <- function(Y, X, conf_x_new_mat, y_grid, theta, sigma, alpha) {
  num_cores <- detectCores() - 1  # Use all but one core
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Export required functions and variables
  clusterExport(cl, varlist = c("posterior_predictive.v2", "Y", "X", "y_grid", "theta", "sigma", "alpha"), envir = environment())
  
  conf_range <- foreach(i = 1:nrow(conf_x_new_mat), .combine = rbind, .packages = c("foreach")) %dopar% {
    x_new_i <- conf_x_new_mat[i, , drop = FALSE]  # Extract the current `x_new`
    
    # Compute conformal prediction set directly within the same loop
    accepted_y <- foreach(l = 1:length(y_grid), .combine = c) %do% {
      # Compute conformity scores on dataset
      sig_1_to_n <- posterior_predictive.v2(Y = Y, X = X, y = y_grid[l], x_new = x_new_i, theta = theta, sigma = sigma)
      
      # Conformity score on test point
      sig_n_plus_one <- posterior_predictive.v2(Y = y_grid[l], X = x_new_i, y = y_grid[l], x_new = x_new_i, theta = theta, sigma = sigma)
      
      # Adjusted quantile calculation
      n <- length(Y)
      pi <- (length(which(sig_1_to_n < sig_n_plus_one)) + 1)
      
      # Only return y_grid[l] if it meets the condition
      if (pi > (floor(alpha*(n+1))+1)) {
        y_grid[l]
      } else {
        NULL
      }
    }
    
    # Compute lower and upper bounds
    bounds <- range(accepted_y)
    return(bounds)
  }
  
  stopCluster(cl)
  
  return(conf_range)  # Matrix where each row is [lower, upper] for an `x_new`
}


#computes posterior predictive via MCMC, AOI importance sampling
posterior_predictive.v2 <- function(Y, X, y, x_new, theta, sigma) {
  # Compute weights(log scale)
  log_w <- apply(theta, 2, function(th) dnorm(y, mean = x_new %*% th, sd = sigma, log=TRUE))
  max_log_w <- max(log_w)
  log_sum_exp <- max_log_w + log(sum(exp(log_w - max_log_w)))
  normalized_weights <- exp(log_w - log_sum_exp)
  
  # Get dimensions
  n_eval <- length(Y)
  n_samples <- ncol(theta)
  
  # Calculate all means for all evaluation points and samples
  all_means <- X %*% theta
  
  # Create an array that repeats Y for each sample
  # Create a n_eval × n_samples matrix where each row is filled with one Y value
  Y_mat <- matrix(Y, nrow=n_eval, ncol=n_samples)
  
  # Calculate all densities at once
  all_densities <- dnorm(Y_mat, mean=all_means, sd=sigma)
  
  # Apply weights and sum rows
  result <- rowSums(all_densities * matrix(normalized_weights, 
                                           nrow=n_eval, 
                                           ncol=n_samples, 
                                           byrow=TRUE))
  
  return(result)
}





#function to get upper and lower normal quantile
credible.bounds <- function(x_new, post.params, alpha)
{
  lower<-qnorm(alpha / 2, mean = x_new%*%post.params$theta_mu, sd = sqrt(sigma^2+x_new%*%post.params$theta_var_mat%*%t(x_new)))
  upper<-qnorm(1-alpha / 2, mean = x_new%*%post.params$theta_mu, sd = sqrt(sigma^2+x_new%*%post.params$theta_var_mat%*%t(x_new)))
  return(c(lower,upper))
}
#compute upper and lower bounds for an x=x-grid
wrapped.credible.bounds <- function(x, post.params, alpha){
# Create a fine grid of x values to approximate the credible bands
cred_x_grid <- seq(min(x), max(x), length.out = 100)
#create design matrix for x_new
cred_x_new_mat<-matrix(c(rep(1,100),cred_x_grid), nrow=100)
cred.bounds <-sapply(1:nrow(cred_x_new_mat), function(i) credible.bounds(cred_x_new_mat[i, , drop=F], post.params, alpha))
return(cred.bounds)
}                                                                       

plot.predicition.ints <- function(features, resp ,conformal.range, credible.range, case="Null", alsamples, besamples){
# Extract lower and upper bounds
conf.lower <- conformal.range[, 1]
conf.upper <-conformal.range[, 2]
cred.l<-credible.range[1,]
cred.u<-credible.range[2,]
# Plot observed data
plot(features, resp, pch = 16, col = "black", xlab = "X", ylab = "Y",
     main = case)
x_grid<-seq(min(features), max(features), length.out = 100)
# Fill the confidence region (blue)
polygon(c(x_grid, rev(x_grid)), c(conf.lower, rev(conf.upper)), col = rgb(0,0,1,0.3), border = NA)

# Fill the region between red dotted lines
polygon(c(x_grid, rev(x_grid)), c(cred.l, rev(cred.u)), 
        col = rgb(1, 0, 0, 0.2), border = NA)

# Plot smooth upper and lower bounds
lines(x_grid, conf.lower, col = cbPalette[5], lwd = 2)
lines(x_grid, conf.upper, col = cbPalette[5], lwd = 2)
lines(x_grid, cred.l, col = cbPalette[7], lwd = 2, lty = 2)  # Lower bound line
lines(x_grid, cred.u, col = cbPalette[7], lwd = 2, lty = 2)  # Upper bound line

# Add the least squares line
# Add the least squares line
abline(a=alsamples$alpha_mean, b=besamples$beta_mean, col =cbPalette[7], lwd = 2)
abline(a = 2, b = 3, lwd = 2)
# Simplified legend with just the two bands
legend("topleft", 
       legend = c("Conformal Band", "Credible Band"),
       lty = c(1, 2),
       lwd = c(2, 2),
       col = c(cbPalette[5], cbPalette[7]),
       fill = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.2)),
       border = NA,
       bg = "white")
}

#finds the average interval lengths for conformal and credible intervals
mean.interval.length <- function(conformal.range,credible.range){
  # Extract lower and upper bounds
  conf.lower <- conformal.range[, 1]
  conf.upper <-conformal.range[, 2]
  cred.l<-credible.range[1,]
  cred.u<-credible.range[2,]
  a<-mean(abs(conf.lower-conf.upper))
  b<-mean(abs(cred.l-cred.u))
  list(average.conf.len=a, average.cred.len=b)
  
}
 

coverage.proportion <- function(x, Y, grid_x, lower_bounds, upper_bounds) {
  # Ensure grid vectors have the same length
  if (length(grid_x) != length(lower_bounds) || length(grid_x) != length(upper_bounds)) {
    stop("grid_x, lower_bounds, and upper_bounds must all have the same length")
  }
  
  # Ensure data vectors have the same length
  if (length(x) != length(Y)) {
    stop("x and Y must have the same length")
  }
  
  # Create interpolation functions for bounds
  lower_interp <- approxfun(grid_x, lower_bounds, rule = 2)
  upper_interp <- approxfun(grid_x, upper_bounds, rule = 2)
  
  # Get interpolated bounds at each data point
  lower_at_x <- lower_interp(x)
  upper_at_x <- upper_interp(x)
  
  # Check which points fall within their interpolated bounds
  inside <- (Y >= lower_at_x) & (Y <= upper_at_x)
  
  list(coverage = mean(inside), 
       threeSE = 3 * sqrt((mean(inside) * (1 - mean(inside))) / length(Y)))
}
#gets max(across parameters) standard error from mcmc samples
gibbs.SE <- function(samples_matrix) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Please install the 'coda' package to compute ESS: install.packages('coda')")
  }
  
  # Check: must be matrix or data frame with >=2 columns
  if (!is.matrix(samples_matrix) && !is.data.frame(samples_matrix)) {
    stop("samples_matrix must be a matrix or data frame with one column per parameter.")
  }
  
  # Avoid variable shadowing — don't use 'samples' inside loop
  se_list <- apply(samples_matrix, 2, function(param_samples) {
    mcmc_obj <- coda::as.mcmc(param_samples)
    param_sd <- sd(param_samples)
    ess <- coda::effectiveSize(mcmc_obj)
    param_sd / sqrt(ess)
  })
  
  max_se <- max(se_list)
  
  return(list(
    max_se = max_se,
    all_se = se_list
  ))
}


# Function to calculate ESS (Effective Sample Size) from MCMC chains
compute_min_ess <- function(theta) {
  # theta should be a matrix with parameters in rows and samples in columns
  
  if (requireNamespace("coda", quietly = TRUE)) {
    library(coda)
    # Convert to mcmc object (coda expects samples in rows, parameters in columns)
    theta_mcmc <- as.mcmc(t(theta))
    # Calculate ESS for each parameter
    ess_values <- effectiveSize(theta_mcmc)
    # Return the minimum ESS
    return(min(ess_values))
    
  } }

# Function to compute and plot ESS with conformal bounds
#usual inputs, select up to 4 alpha values
#mcmc_ess is for scaling plot by min(ESS)/T
#linelen is a plotting parameter to select length to wich conformla bound will be extended vertically
plot_ess_with_bounds <- function(Y, X, y_grid, x_new, theta, sigma, 
                                 mcmc_ess = NULL, 
                                 alpha_values = c(0.1, 0.2),
                                 title, linelen) {
  
  require(doParallel)
  require(foreach)
  
  n_grid <- length(y_grid)
  ess_values <- numeric(n_grid)
  n_samples <- ncol(theta)
  
  # Compute ESS for each value of y in the grid
  for (i in 1:n_grid) {
    y <- y_grid[i]
    
    # Compute weights (still using log scale for numerical stability) 
    log_w <- apply(theta, 2, function(th) dnorm(y, mean = x_new %*% th, sd = sigma, log=TRUE))
    max_log_w <- max(log_w)
    log_sum_exp <- max_log_w + log(sum(exp(log_w - max_log_w)))
    normalized_weights <- exp(log_w - log_sum_exp)
    
    # Calculate ESS using formula 1/sum(w_i^2) for normalized weights
    ess_values[i] <- 1 / sum(normalized_weights^2)
  }
  
  # Scale ESS values 
  if (!is.null(mcmc_ess)) {
    ess_values <- ess_values *(mcmc_ess/  n_samples)
  }
  
  # Create plot
  plot(y_grid, ess_values, type="l", lwd=2,
       xlab="y", ylab="Effective Sample Size",
       main=title)
  
  # Format x_new for conformal_bounds_single_parallel
  conf_x_new_mat <- x_new
  
  # Define colors for different alpha values
  alpha_colors <- c("red", "blue", "green", "purple", "orange")[1:length(alpha_values)]
  
  # Compute and plot bounds for each alpha
  for (i in 1:length(alpha_values)) {
    alpha <- alpha_values[i]
    
    # Compute conformal bounds
    conf_bounds <- conformal_bounds_single_parallel(
      Y = Y, X = X, 
      conf_x_new_mat = conf_x_new_mat, 
      y_grid = y_grid, 
      theta = theta, 
      sigma = sigma, 
      alpha = alpha
    )
    
    # Extract lower and upper bounds
    lower_bound <- conf_bounds[1]
    upper_bound <- conf_bounds[2]
    # Find y-positions for segments near the ESS curve
    # Find the closest y_grid points to our bounds
    lower_idx <- which.min(abs(y_grid - lower_bound))
    upper_idx <- which.min(abs(y_grid - upper_bound))
    
    # Get the ESS values at these points
    lower_ess <- ess_values[lower_idx]
    upper_ess <- ess_values[upper_idx]
    
    # Add localized vertical line segments near the ESS curve
    segments(
      x0 = lower_bound, y0 = lower_ess-linelen,
      x1 = lower_bound, y1 = lower_ess+linelen,  # Extend slightly above the curve
      col = alpha_colors[i], lwd = 2
    )
    segments(
      x0 = upper_bound, y0 = upper_ess-linelen,
      x1 = upper_bound, y1 = upper_ess+linelen,  # Extend slightly above the curve
      col = alpha_colors[i], lwd = 2
    )
    
  }
  
  # Add legend
  legend("topright", 
         legend = sapply(alpha_values, function(a) bquote(alpha == .(a))),
         col = alpha_colors, 
         lty = 1, 
         lwd = 2,
         cex = 0.8)
  
}

