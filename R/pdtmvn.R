#' Function to evaluate the density of a pdTMVN distribution
#' 
#' @param x Vector or matrix of quantiles.  If x is a matrix, each row is taken
#'   to be a quantile.
#' @param mean Mean vector, default is rep(0, nrow(sigma)).
#' @param sigma Covariance matrix, default is diag(length(mean)).
#' @param precision Precision (inverse covariance) matrix.
#' @param lower Vector of lower truncation points, default is
#'   rep(-Inf, length = length(mean)).
#' @param upper Vector of upper truncation points, default is 
#'   rep(Inf, length = length(mean)).
#' @param norm_const NOT CURRENTLY USED.  If many calls will be made to this
#'   function with the same covariance/precision, it may be helpful to
#'   precompute the normalization constant for the distribution and pass that in
#'   so that it is not re-calculated on every function call
#' @param trunc_const If many calls will be made to this function with the same
#'   covariance/precision and lower/upper truncation bounds, it may be helpful
#'   to precompute the truncation constant for the distribution and pass that
#'   in so that it is not re-calculated on every function call
#' @param sigma_continuous Matrix containing the marginal
#'   covariance matrix for the continuous variables.  Providing this saves a
#'   small amount of computation time if precision is provided but sigma isn't.
#' @param conditional_sigma_discrete Matrix containing the conditional covariance
#'   of the discrete variables given the continuous variables.  Providing this
#'   saves a small amount of computation time if sigma is provided but precision
#'   isn't.
#' @param conditional_mean_discrete_offset_multiplier Matrix containing
#'   Sigma_dc Sigma_c^{-1}.  This is used in computing the mean of the underlying
#'   multivariate normal distribution for the discrete variables conditioning on
#'   the continuous variables.  Providing this saves a small amount of
#'   computation time.
#' @param continous_vars Vector containing either column names or column
#'   indices for continous variables.
#' @param discrete_vars Vector containing either column names or column indices
#'   for discrete variables.
#' @param discrete_var_range_functions a list with one entry for each element
#'   of discrete_vars.  Each entry is a character vector of length 2; the first
#'   element is the name of a function that returns a(x) for any real x, and
#'   the second is the name of a function that returns b(x) for any real x.
#' @param log Logical; if TRUE, return value is log(density).
#' @param validate_level Numeric; if 0, no validation is performed (risky!).
#'   If 1, parameters are validated but warnings about checks not performed are
#'   not issued.  If 2, parameters are validated and warnings about checks not
#'   performed are issued.
#' 
#' @return Density of pdTMVN distribution evaluated at x
dpdtmvn <- function(x,
		mean = rep(0, nrow(sigma)),
		sigma = diag(length(mean)),
		precision,
		lower = rep(-Inf, length = length(mean)),
		upper = rep(Inf, length = length(mean)),
		norm_const,
		trunc_const,
		sigma_continuous,
		conditional_sigma_discrete,
		conditional_mean_discrete_offset_multiplier,
		continuous_vars = seq_along(mean),
		discrete_vars = NULL,
		discrete_var_range_functions = sapply(seq_along(discrete_vars), function(dv) {
			c("floor_x_minus_1", "floor")
		}),
		log = FALSE,
		validate_level = TRUE) {
	
	## Convert x to matrix if a vector or data frame was passed in
	if(is.vector(x)) {
		x_names <- names(x)
		dim(x) <- c(1, length(x))
		colnames(x) <- x_names
	} else if(is.data.frame(x)) {
		x <- as.matrix(x)
	}
	
	## Validate parameters
	if(validate_level > 0) {
		validated_params <- validate_params_pdtmvn(x = x,
			mean = mean,
			sigma = sigma,
			precision = precision,
			lower = lower,
			upper = upper,
			norm_const = norm_const,
			trunc_const = trunc_const,
			sigma_continuous = sigma_continuous,
			conditional_sigma_discrete = conditional_sigma_discrete,
			conditional_mean_discrete_offset_multiplier = 
				conditional_mean_discrete_offset_multiplier,
			continuous_vars = continuous_vars,
			discrete_vars = discrete_vars,
			validate_level = validate_level)
		
		## validate_params_pdtmvn may update some values that were passed in.
		## Assign those updated values to the corresponding variables in the
		## current environment.
		for(var_name in names(validated_params)) {
			assign(var_name, validated_params[[var_name]])
		}
		
		rm("validated_params")
	}
	
	## The outline for the rest of the code here is closely based on that in the
	## tmvtnorm package, but we
	##  - do the integration for partial discretization
	##  - allow for precomputation of trunc_const
	##    (and possibly norm_const in the future)
	
	## determine whether each observation is in the support of the distribution
	## 2 ways an x vector could fail to be in the support:
	##  - a discrete or continuous covariate falls above the upper truncation bound
	##      or below the lower truncation bound for that covariate
	##  - a discrete covariate is not at one of the points where the discrete
	##      distribution for that covariate places positive mass.
	
	## truncation limits
	in_truncation_support <- which(apply(x, 1, function(x_row) {
		all(x[i, ] >= lower & x[i, ] <= upper & !any(is.infinite(x)))
	}))
	
	## discrete distribution domain
	if(length(discrete_vars) > 0) {
		## if there were a lot of x's outside of truncation support, we could save
		## time here -- but for now, I expect such values to be rare.
		
		## the vector b_x for each row of x
		b_x_discrete <- sapply(seq_along(discrete_vars), function(discrete_var_ind) {
			do.call(discrete_var_range_functions[[discrete_var_ind]][[2]],
				list(x=x[, discrete_vars[discrete_var_ind]])
			)
		})
		
		## logical vector of length nrow(x) with whether all entries corresponding to
		## discrete variables in row i of x agree with all entries of row i of
		## b_x_discrete
		in_discrete_dist_domain <- sapply(seq_len(nrow(x)),
			function(x_row_ind) {
				isTRUE(all.equal(b_x_discrete[x_row_ind], x[x_row_ind, discrete_vars]))
			})
	} else {
		in_discrete_dist_domain <- rep(TRUE, nrow(x))
	}
	
	## in_support are now row indices for observations in support
	in_support <- which(in_truncation_support & in_discrete_dist_domain)
	
	## compute result -- computations are on log scale
	log_result <- rep(-Inf, nrow(x))
	
	## Compute contribution from observations of continuous variables
	if(length(continuous_vars) > 0) {
		log_result[in_support] <- dmvnorm(x[in_support, continuous_vars, drop=FALSE],
			mean = mean[continuous_vars],
			sigma = sigma_continuous_vars,
			log = TRUE)
	}
	
	## Compute contribution from observations of discrete variables
	if(length(discrete_vars) > 0) {
		## get conditional means of discrete vars given continuous vars
		cond_means <- matrix(rep(mean[discrete_vars], nrow(x)),
			nrow=nrow(x),
			byrow=TRUE)
		
		if(length(continuous_vars) > 0) {
			## compute x_ci - mu_c, stacked in a matrix
			x_c_differences <- sweep(x[in_support, continuous_vars, drop=FALSE],
				2,
				mean[continuous_vars],
				`-`)
			
			## adjust -- this step gives us a matrix where row i is the transpose of
			## mu_d + Sigma_{dc} Sigma_{c}^{-1} (x_ci - mu_c)
			cond_means <- cond_means +
				x_c_differences %*% t(conditional_mean_discrete_offset_multiplier)
		}
		
		## get lower bounds on integration for each discrete observation
		a_x_discrete <- sapply(seq_along(discrete_vars), function(discrete_var_ind) {
			do.call(discrete_var_range_functions[[discrete_var_ind]][[1]],
							list(x=x[in_support, discrete_vars[discrete_var_ind]])
			)
		})
		
		## we computed b_x earlier for all x -- subset here for those x in support
		b_x_discrete <- b_x_discrete[in_support, , drop=FALSE]
		
		if(identical(length(discrete_vars), 1L)) {
			p_lt_b <- pnorm(b_x_discrete[in_support, discrete_vars],
				mean = as.vector(cond_means),
				sd = sqrt(as.vector(conditional_sigma_discrete)))
			
			p_lt_a <- pnorm(a_x_discrete[in_support, discrete_vars],
				mean = as.vector(cond_means),
				sd = sqrt(as.vector(conditional_sigma_discrete)))
			
			log_result[in_support] <- log_result[in_support] +
				logspace_sub(p_lt_b, p_lt_a)
		} else {
			stop("dpdtmvn does not currently support discrete_vars with length > 1")
			## pmvnorm does not support a log=TRUE option; I suspect that we need that
			## for the probabilities to be non-zero.
			log_result[in_support] <- log_result[in_support] +
				apply(matrix(seq_along(in_support)), 1, function(support_row_ind) {
					pmvnorm(lower=a_x_discrete[support_row_ind, ],
									upper=b_x_discrete[support_row_ind, ],
									mean=mean[discrete_vars],
									sigma=conditional_sigma_discrete)
				})
		}
	}
	
	## Adjust for truncation by subtracting trunc_const.
	log_result[in_support] <- log_result[in_support] - trunc_const
		
	## Return
	if(log) {
		return(log_result)
	} else {
		return(exp(log_result))
	}
}
	
#' Validate the parameters of a call to dpdtmvn
#' 
#' @param x matrix of quantiles.  If x is a matrix, each row is taken to be a
#'   quantile.
#' @param mean Mean vector, default is rep(0, nrow(sigma)).
#' @param sigma Covariance matrix, default is diag(length(mean)).
#' @param precision Precision (inverse covariance) matrix.
#' @param lower Vector of lower truncation points, default is
#'   rep(-Inf, length = length(mean)).
#' @param upper Vector of upper truncation points, default is 
#'   rep(Inf, length = length(mean)).
#' @param norm_const If many calls will be made to this function with the same
#'   covariance/precision, it may be helpful to precompute the normalization
#'   constant for the distribution and pass that in so that it is not
#'   re-calculated on every function call
#' @param trunc_const If many calls will be made to this function with the same
#'   covariance/precision and lower/upper truncation bounds, it may be helpful
#'   to precompute the truncation constant for the distribution and pass that
#'   in so that it is not re-calculated on every function call
#' @param sigma_continuous Matrix containing the marginal
#'   covariance matrix for the continuous variables.  Providing this saves a
#'   small amount of computation time if precision is provided but sigma isn't.
#' @param conditional_sigma_discrete Matrix containing the conditional covariance
#'   of the discrete variables given the continuous variables.  Providing this
#'   saves a small amount of computation time if sigma is provided but precision
#'   isn't.
#' @param conditional_mean_discrete_offset_multiplier Matrix containing
#'   Sigma_dc Sigma_c^{-1}.  This is used in computing the mean of the underlying
#'   multivariate normal distribution for the discrete variables conditioning on
#'   the continuous variables.  Providing this saves a small amount of
#'   computation time.
#' @param continous_vars Vector containing either column names or column indices for
#'   continous variables.
#' @param discrete_vars Vector containing either column names or column indices for
#'   discrete variables.
#' @param validate_level Numeric; if <= 1, validation is performed, but warnings
#'   about things that were not checked are not issued.  If 2, validation is
#'   performed and warnings about things that were not checked are issued.
#' 
#' @return Named list with parameter values after validation.
validate_params_pdtmvn <- function(x,
	mean,
	sigma,
	precision,
	lower,
	upper,
	norm_const,
	trunc_const,
	sigma_continuous,
	conditional_sigma_discrete,
	conditional_mean_discrete_offset_multiplier,
	continuous_vars,
	discrete_vars,
	validate_level) {
	
	if(validate_level > 1) {
		warning("parameter validation does not yet check whether your truncation limits make sense relative to your functions for computing the range of you discrete random variables.")
		warning("parameter validation does not yet check whether your discrete_var_range_functions are in the same order as the discrete columns of x")
	}
	
	## list of results
	validated_params <- list()
	
	
	## Validate x
	if(!missing(x)) {
		if(!is.numeric(x)) {
			stop("x must be a numeric matrix, vector, or data frame")
		}
		validated_params$x <- x
	}
	
	
	## Validate mean
	if(!missing(mean)) {
		if(!is.numeric(mean) || !identical(length(mean), ncol(x))) {
			stop("mean parameter must be a numeric vector with the same length as the number of covariates in x.")
		}
		validated_params$mean <- mean
	}
	
	
	## Validate precision/covariance
	if(!missing(sigma)) {
		if(!is.numeric(sigma) || !is.matrix(sigma) ||
			 	!(identical(nrow(sigma), ncol(x)) &&
			  	identical(ncol(sigma), ncol(x)))) {
			stop("sigma parameter must be a numeric J by J matrix where J is the number of covariates in x.")
		}
		if(validate_level > 1) {
			warning("did not check that sigma is a p.d. symmetric matrix")
		}
		validated_params$sigma <- sigma
	}
	if(!missing(precision)) {
		if(!is.numeric(precision) || !is.matrix(precision) ||
			 	!(identical(nrow(precision), ncol(x)) &&
			 			identical(ncol(precision), ncol(x)))) {
			stop("precision parameter must be a numeric J by J matrix where J is the number of covariates in x.")
		}
		validated_params$precision <- precision
		if(validate_level > 1) {
			warning("did not check that precision is a p.d. symmetric matrix")
			if(!missing(sigma)) {
				warning("did not check that precision and sigma agree")
			}
		}
	}
	if(missing(sigma) && missing(precision)) {
		validated_params$sigma <- sigma ## fill with default value
	}
	
	
  ## Validate lower and upper bounds
	if(!missing(lower)) {
		if(!is.numeric(lower) || !identical(length(lower), ncol(x))) {
			stop("lower parameter must be a numeric vector with the same length as the number of covariates in x.")
		}
		validated_params$lower <- lower
	}
	if(!missing(upper)) {
		if(!is.numeric(upper) || !identical(length(upper), ncol(x))) {
			stop("upper parameter must be a numeric vector with the same length as the number of covariates in x.")
		}
		validated_params$upper <- upper
	}
	
	
	## Validate norm_const and trunc_const; compute them if they are missing
	if(!missing(norm_const)) {
		if(!is.numeric(norm_const) || !identical(length(norm_const), 1L)) {
			stop("norm_const parameter must be numeric with length 1.")
		}
	} else {
		## The value of norm_const is not currently used -- implement this calculation
		## if it is used at some point in the future.
#		norm_const <- compute_norm_const_pdtmvn(sigma, precision)
		norm_const <- NA
	}
	validated_params$norm_const <- norm_const
	
	if(!missing(trunc_const)) {
		if(!is.numeric(trunc_const) || !identical(length(trunc_const), 1L)) {
			stop("trunc_const parameter must be numeric with length 1.")
		}
	} else {
		trunc_const <- compute_trunc_const_pdtmvn(mean = mean,
			sigma = sigma,
			precision = precision,
			lower = lower,
			upper = upper)
	}
	validated_params$trunc_const <- trunc_const
	
	
	## Validate continous and discrete_vars.  We do the following:
	##  1) convert character vectors to numeric.
	##  2) (a) If both were specified, make sure they form a partition of the
	##     column indices of x
	##     (b) If only one was specified, set the other to its complement
	if(!missing(continuous_vars)) {
		if(is.character(continous_vars)) {
			continuous_vars <- which(colnames(x) %in% continuous_vars)
		} else if(is.numeric(continuous_vars)) {
			continuous_vars <- as.integer(continuous_vars)
		} else {
			stop("continuous_vars parameter must be a character or numeric vector.")
		}
		
		## Initial check of data type for discrete_vars if it was provided.
		if(!missing(discrete_vars)) {
			if(is.character(discrete_vars)) {
				discrete_vars <- which(colnames(x) %in% discrete_vars)
			} else if(is.numeric(discrete_vars)) {
				discrete_vars <- as.integer(discrete_vars)
			} else {
				stop("continuous_vars parameter must be a character or numeric vector.")
			}
			
			## Both continuous_vars and discrete_vars were specified;
			## ensure that they specify a partition (disjoint cover)
			## of the columns of x
			all_specified_vars <- c(continuous_vars, discrete_vars)
			if(!all(seq_len(ncol(x)) %in% all_specified_vars) ||
				 	!identical(length(all_specified_vars), ncol(x))) {
				stop("discrete_vars and continuous_vars must include all columns of x without any overlap.")
			}
		} else {
			discrete_vars <- which(!(seq_len(ncol(x)) %in% continuous_vars))
		}
	} else {
		## continuous_vars is missing
		## If discrete_vars is non-missing, override the default value of
		## continuous_vars specified in the function header to set continuous_vars
		## equal to the complement of discrete_vars.
		## Otherwise, if discrete_vars is also missing, do nothing.
		
		## Initial check of data type for discrete_vars if it was provided.
		if(!missing(discrete_vars)) {
			if(is.character(discrete_vars)) {
				discrete_vars <- which(colnames(x) %in% discrete_vars)
			} else if(is.numeric(discrete_vars)) {
				discrete_vars <- as.integer(discrete_vars)
			} else {
				stop("continuous_vars parameter must be a character or numeric vector.")
			}
			
			## Set continuous_vars to the complement of discrete_vars.
			continuous_vars <- which(!(seq_len(ncol(x)) %in% discrete_vars))
		}
	}

	if(length(discrete_vars) > 1) {
		stop("dpdtmvn does not currently support discrete_vars with length > 1.")
	}

	validated_params$continuous_vars <- continuous_vars
	validated_params$discrete_vars <- discrete_vars
	
	
	## fill in values for sigma_continuous, conditional_sigma_discrete, and
	## conditional_mean_discrete_offset_multiplier if they are missing
	## if they are provided, these parameter values are not validated
	if(length(continuous_vars) > 0) {
		## Get sigma_continous
		if(missing(sigma_continuous)) {
			if("sigma" %in% names(validated_params)) {
				## sigma_continuous is a subset of sigma
				validated_params$sigma_continuous <- 
					sigma[continuous_vars, continuous_vars, drop=FALSE]
			} else {
				## sigma_continuous can be calculated from precision
				if(identical(length(discrete_vars), 0L)) {
					validated_params$sigma_continuous <- solve(precision)
				} else {
					validated_params$sigma_continuous <- calc_Schur_complement(precision,
						continuous_vars)
				}
			}
		} else {
			validated_params$sigma_continuous <- sigma_continuous
			if(validate_level > 1) {
				warning("did not check that provided sigma_continuous agrees with sigma/precision")
			}
		}
		
		## if there are both continuous and discrete variables, get
		## conditional_sigma_discrete and conditional_mean_discrete_offset_multiplier
		if(length(discrete_vars) > 0) {
			## Get conditional_sigma_discrete
			if(missing(conditional_sigma_discrete)) {
				if(!missing(precision)) {
					precision_d <- precision[discrete_vars, discrete_vars, drop=FALSE]
					validated_params$conditional_sigma_discrete <- solve(precision_d)
				} else {
					validated_params$conditional_sigma_discrete <- 
						calc_Schur_complement(validated_params$sigma, discrete_vars)
				}
			} else {
				validated_params$conditional_sigma_discrete <- conditional_sigma_discrete
				if(validate_level > 1) {
					warning("did not check that provided conditional_sigma_discrete agrees with sigma/precision")
				}
			}
			
			## Get conditional_mean_discrete_offset_multiplier
			if(missing(conditional_mean_discrete_offset_multiplier)) {
				if(!missing(sigma)) {
					sigma_dc <- sigma[discrete_vars, continuous_vars, drop=FALSE]
					
					validated_params$conditional_mean_discrete_offset_multiplier <-
						sigma_dc * solve(validated_params$sigma_continuous)
				} else {
					precision_d <- precision[discrete_vars, discrete_vars, drop=FALSE]
					precision_dc <- precision[discrete_vars, continuous_vars, drop=FALSE]

					validated_params$conditional_mean_discrete_offset_multiplier <-
						-1 * solve(precision_d) %*% precision_dc
				}
			} else {
				validated_params$conditional_mean_discrete_offset_multiplier <-
					conditional_mean_discrete_offset_multiplier
				if(validate_level > 1) {
					warning("did not check that provided conditional_mean_discrete_offset_multiplier agrees with sigma/precision")
				}
			}
		}
	} else if(length(discrete_vars) > 0) {
		## only discrete variables -- we only need conditional_sigma_discrete,
		## not sigma_continuous or conditional_mean_discrete_offset_multiplier
		if(missing(conditional_sigma_discrete)) {
			if(!missing(sigma)) {
				validated_params$conditional_sigma_discrete <- sigma
			} else {
				validated_params$conditional_sigma_discrete <- solve(precision)
			}
		} else {
			validated_params$conditional_sigma_discrete <- conditional_sigma_discrete
			if(validate_level > 1) {
				warning("did not check that provided conditional_sigma_discrete agrees with sigma/precision")
			}
		}
	}
	
	
	## Return
	return(list(x = x,
	 mean = mean,
	 sigma = sigma,
	 precision = precision,
	 lower = lower,
	 upper = upper,
	 norm_const = norm_const,
	 trunc_const = trunc_const,
	 continuous_vars = continuous_vars,
	 discrete_vars = discrete_vars))
}


compute_trunc_const_pdtmvn <- function(mean, sigma, precision, lower, upper) {
	if(missing(sigma)) {
		sigma <- solve(precision)
	}
	
	## Ideally, we would like pmvnorm to support a log=TRUE argument to give us
	## more numerical precision in high dimensional cases.  Here we compute the
	## probability and then take its logarithm.
	exp_trunc_const <- pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
	
	return(log(exp_trunc_const))
}


floor_x_minus_1 <- function(x) {
	return(floor(x) - 1)
}

calc_Schur_complement <- function(X, inds) {
	X_1 <- precision[continuous_vars, continuous_vars, drop=FALSE]
	X_12 <- precision[continuous_vars, discrete_vars, drop=FALSE]
	X_2 <- precision[discrete_vars, discrete_vars, drop=FALSE]
	
	return(solve(X_1 - X_12 %*% solve(X_2) %*% t(X_12)))
}
