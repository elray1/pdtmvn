## Miscellaneous functions used in evaluating and sampling from the pdtmvn
## distribution.
## 
## in_pdtmvn_support
## calc_in_truncation_support
## validate_params_pdtmvn
## compute_sigma_subcomponents
## get_conditional_mvn_intermediate_params
## get_conditional_mvn_params
## get_conditional_mvn_mean_from_intermediate_params
## compute_trunc_const_pdtmvn
## equals_integer
## floor_x_minus_1
## calc_Schur_complement

#' Function to determine whether observations are within the support of a
#' pdTMVN distribution
#' 
#' @param x Vector or matrix of quantiles.  If x is a matrix, each row is taken
#'   to be a quantile.
#' @param lower Vector of lower truncation points, default is
#'   rep(-Inf, length = length(mean)).
#' @param upper Vector of upper truncation points, default is 
#'   rep(Inf, length = length(mean)).
#' @param continuous_vars Vector containing either column names or column
#'   indices for continuous variables.
#' @param discrete_vars Vector containing either column names or column indices
#'   for discrete variables.
#' @param discrete_var_range_fns a list with one entry for each element
#'   of discrete_vars.  Each entry is a named list of length 3; the element
#'   named "a" is a character string with the name of a function that returns
#'   a(x) for any real x, the element named "b" is a character string with the
#'   name of a function that returns b(x) for any real x, and the element named
#'   "in_range" is a character string with the name of a function that returns a
#'   logical, TRUE if x is in the support of the corresponding discrete variable
#'   and FALSE otherwise.
#' 
#' @return Integer vector with indices of rows of x that are in the support of
#'   the pdtmvn distribution
in_pdtmvn_support <- function(x,
        lower = rep(-Inf, length = ncol(x)),
        upper = rep(Inf, length = ncol(x)),
        continuous_vars,
        discrete_vars,
        discrete_var_range_fns) {
    ## determine whether each observation is in the support of the distribution
    ## 2 ways an x vector could fail to be in the support:
    ##  - a discrete or continuous covariate falls above the upper truncation bound
    ##      or below the lower truncation bound for that covariate
    ##  - a discrete covariate is not at one of the points where the discrete
    ##      distribution for that covariate places positive mass.
    
    ## truncation limits
	in_truncation_support <- calc_in_truncation_support(x, lower, upper)
	
	## discrete distribution domain
	if(length(discrete_vars) > 0) {
		## if there were a lot of x's outside of truncation support, we could save
		## time here -- but for now, I expect such values to be rare.
		
		## the vector b_x for each row of x
		b_x_discrete <- plyr::laply(seq_along(discrete_vars), function(discrete_var_ind) {
			do.call(discrete_var_range_fns[[discrete_var_ind]][["b"]],
				list(x=x[, discrete_vars[discrete_var_ind]])
			)
		})
		if(identical(length(discrete_vars), 1L)) {
			b_x_discrete <- matrix(b_x_discrete)
		} else {
			b_x_discrete <- t(b_x_discrete)
		}
		
		## logical vector of length nrow(x) with whether all entries corresponding to
		## discrete variables in row i of x are in their domains
		in_discrete_dist_domain <- plyr::laply(seq_along(discrete_vars), function(discrete_var_ind) {
		    do.call(discrete_var_range_fns[[discrete_var_ind]][["in_range"]],
		        list(x=x[, discrete_vars[discrete_var_ind]])
		    )
		})
		if(length(discrete_vars) > 1L) {
		    in_discrete_dist_domain <- apply(t(in_discrete_dist_domain), 1, all)
		}
	} else {
		in_discrete_dist_domain <- rep(TRUE, nrow(x))
	}
	
	## in_support are now row indices for observations in support
	in_support <- which(in_truncation_support & in_discrete_dist_domain)
    
    return(in_support)
}

#' Function to determine whether observations are within the truncation support of a
#' pdTMVN distribution
#' 
#' @param x Matrix of quantiles.  If x is a matrix, each row is taken
#'   to be a quantile.
#' @param lower Vector of lower truncation points
#' @param upper Vector of upper truncation points
#' 
#' @return Logical vector with length = number of rows of x, where entry i is TRUE if row i
#'   is between the bounds specified by lower and upper and FALSE otherwise
calc_in_truncation_support <- function(x, lower, upper) {
	storage.mode(x) <- "double"
	storage.mode(lower) <- "double"
	storage.mode(upper) <- "double"
	
	return(as.logical(.Call("in_truncation_support_C", x, lower, upper)))
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
#' @param continuous_vars Vector containing either column names or column indices for
#'   continuous variables.
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
		warning("parameter validation does not yet check whether your discrete_var_range_fns are in the same order as the discrete columns of x")
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
	} else {
		validated_params$mean <- rep(0, nrow(sigma))
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
		validated_params$sigma <- rep(0, nrow(sigma)) ## fill with default value
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
	
	
	## Validate continuous and discrete_vars.  We do the following:
	##  1) convert character vectors to numeric.
	##  2) (a) If both were specified, make sure they form a partition of the
	##     column indices of x
	##     (b) If only one was specified, set the other to its complement
	if(!missing(continuous_vars) && !is.null(continuous_vars) && length(continuous_vars) > 0) {
		if(is.character(continuous_vars)) {
			continuous_vars <- which(colnames(x) %in% continuous_vars)
		} else if(is.numeric(continuous_vars)) {
			continuous_vars <- as.integer(continuous_vars)
		} else {
			stop("continuous_vars parameter must be a character or numeric vector.")
		}
		
		## Initial check of data type for discrete_vars if it was provided.
		if(!missing(discrete_vars) && !is.null(discrete_vars) && length(discrete_vars) > 0) {
			if(is.character(discrete_vars)) {
				discrete_vars <- which(colnames(x) %in% discrete_vars)
			} else if(is.numeric(discrete_vars)) {
				discrete_vars <- as.integer(discrete_vars)
			} else {
				stop("discrete_vars parameter must be a character or numeric vector.")
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
		## Otherwise, if discrete_vars is also missing, set all vars to continuous.
		
		## Initial check of data type for discrete_vars if it was provided.
		if(!missing(discrete_vars) && !is.null(discrete_vars) && length(discrete_vars) > 0) {
			if(is.character(discrete_vars)) {
				discrete_vars <- which(colnames(x) %in% discrete_vars)
			} else if(is.numeric(discrete_vars)) {
				discrete_vars <- as.integer(discrete_vars)
			} else {
				stop("continuous_vars parameter must be a character or numeric vector.")
			}
			
			## Set continuous_vars to the complement of discrete_vars.
			continuous_vars <- which(!(seq_len(ncol(x)) %in% discrete_vars))
		} else {
			continuous_vars <- seq_len(ncol(x))
			discrete_vars <- integer(0)
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
    validated_params <- c(validated_params,
        compute_sigma_subcomponents(sigma = validated_params$sigma,
            precision = validated_params$precision,
            continuous_vars = continuous_vars,
            discrete_vars = discrete_vars,
            sigma_continuous = sigma_continuous,
            conditional_sigma_discrete = conditional_sigma_discrete,
            conditional_mean_discrete_offset_multiplier = 
                conditional_mean_discrete_offset_multiplier,
            validate_level = validate_level))
	
	## Return
	return(validated_params)
}

#' Compute sigma_continuous, conditional_sigma_discrete, and
#' conditional_mean_discrete_offset_multiplier from sigma and/or precision.
#' These quantities are used in computation of the pdtmvn density.
#' 
#' @param sigma the covariance matrix parameter for the pdtmvn distribution
#' @param precision the precision (inverse covariance) matrix parameter for the
#'     pdtmvn distribution
#' @param continuous_vars integer vector with indices of sigma/precision that
#'     correspond to continuous variables in the pdtmvn distribution.
#' @param discrete_vars integer vector with indices of sigma/precision that
#'     correspond to discrete variables in the pdtmvn distribution.
#' @param sigma_continuous (optional) this matrix may be provided.  If it is,
#'     the provided value will be included in the return value; the value will
#'     not be validated.
#' @param conditional_sigma_discrete (optional) this matrix may be provided.  If
#'     it is, the provided value will be included in the return value; the value
#'     will not be validated.
#' @param conditional_mean_discrete_offset_multiplier (optional) this matrix may
#'     be provided.  If it is, the provided value will be included in the return
#'     value; the value will not be validated.
#' @param validate_level Numeric; if <= 1, validation is performed, but warnings
#'   about things that were not checked are not issued.  If 2, validation is
#'   performed and warnings about things that were not checked are issued.
#'     
#' @return a named list with three components: sigma_continuous,
#'     conditional_sigma_discrete, and
#'     conditional_mean_discrete_offset_multiplier
compute_sigma_subcomponents <- function(sigma = NULL,
    precision = NULL,
    continuous_vars,
    discrete_vars,
    sigma_continuous,
    conditional_sigma_discrete,
    conditional_mean_discrete_offset_multiplier,
    validate_level) {
    if(is.null(sigma) && is.null(precision)) {
        stop("At least one of sigma and precision must be provided.")
    }
    
    sigma_subcomponents <- list()
    
    if(length(continuous_vars) > 0) {
        ## Get sigma_continuous
        if(missing(sigma_continuous)) {
            if(!is.null(sigma)) {
                ## sigma_continuous is a subset of sigma
                sigma_subcomponents$sigma_continuous <- 
                    sigma[continuous_vars, continuous_vars, drop=FALSE]
            } else {
                ## sigma_continuous can be calculated from precision
                if(identical(length(discrete_vars), 0L)) {
                    sigma_subcomponents$sigma_continuous <- solve(precision)
                } else {
                    sigma_subcomponents$sigma_continuous <- 
                        solve(calc_Schur_complement(precision, continuous_vars))
                }
            }
        } else {
            sigma_subcomponents$sigma_continuous <- sigma_continuous
            if(validate_level > 1) {
                warning("did not check that provided sigma_continuous agrees with sigma/precision")
            }
        }
        
        ## if there are both continuous and discrete variables, get
        ## conditional_sigma_discrete and conditional_mean_discrete_offset_multiplier
        if(length(discrete_vars) > 0) {
            temp <- get_conditional_mvn_intermediate_params(sigma = sigma,
                precision = precision,
                fixed_vars = continuous_vars,
                free_vars = discrete_vars,
                conditional_sigma = conditional_sigma_discrete,
                conditional_mean_offset_multiplier = conditional_mean_discrete_offset_multiplier,
                validate_level = validate_level)
            
            sigma_subcomponents$conditional_sigma_discrete <-
                temp$conditional_sigma
            sigma_subcomponents$conditional_mean_discrete_offset_multiplier <-
                temp$conditional_mean_offset_multiplier
        }
    } else if(length(discrete_vars) > 0) {
        ## only discrete variables -- we only need conditional_sigma_discrete,
        ## not sigma_continuous or conditional_mean_discrete_offset_multiplier
        if(missing(conditional_sigma_discrete)) {
            if(!is.null(sigma)) {
                sigma_subcomponents$conditional_sigma_discrete <- sigma
            } else {
                sigma_subcomponents$conditional_sigma_discrete <-
                    solve(precision)
            }
        } else {
            sigma_subcomponents$conditional_sigma_discrete <-
                conditional_sigma_discrete
            if(validate_level > 1) {
                warning("did not check that provided conditional_sigma_discrete agrees with sigma/precision")
            }
        }
    }
    
    return(sigma_subcomponents)
}

#' Compute "intermediate paramaters" for the conditional distribution of a
#' subset of variables in a multivariate normal random vector given the other
#' variables in the vector.  These are the conditional covariance matrix and the
#' matrix that is used as a multiplier to obtain the conditional mean.  As a
#' convenience, either of these values may be supplied, in which case they are
#' returned without any checks of their validity.
#' 
#' @param sigma covariance matrix for the full random vector
#' @param precision precision matrix for the full random vector
#' @param fixed_vars vector of names or column/row indices corresponding to
#'   variables that we are conditioning on
#' @param free_vars vector of names or column/row indices corresponding to
#'   variables whose conditional distribution we want
#' @param conditional_sigma (OPTIONAL) conditional covariance matrix of the
#'   free variables given the fixed variables
#' @param conditional_mean_offset_multiplier (OPTIONAL) matrix used in computing
#'   the conditional mean of the free variables given the fixed variables
#' 
#' @return a named list with two entries: conditional_sigma and
#'   conditional_mean_offset_multiplier
get_conditional_mvn_intermediate_params <- function(sigma,
    precision,
    fixed_vars,
    free_vars,
    conditional_sigma,
    conditional_mean_offset_multiplier,
    validate_level = 1L) {
    
    result <- list()
    
    ## Get conditional_sigma
    if(missing(conditional_sigma) || is.null(conditional_sigma)) {
        if(!missing(precision) && !is.null(precision)) {
            precision_free <- precision[free_vars, free_vars, drop=FALSE]
            result$conditional_sigma <- solve(precision_free)
        } else {
            result$conditional_sigma <- calc_Schur_complement(sigma, free_vars)
        }
    } else {
        result$conditional_sigma <- conditional_sigma
        if(validate_level > 1) {
            warning("did not check that provided conditional_sigma_discrete agrees with sigma/precision")
        }
    }
    
    ## Get conditional_mean_offset_multiplier
    if(length(fixed_vars) > 0) {
        if(missing(conditional_mean_offset_multiplier)) {
            if(!is.null(sigma)) {
                sigma_fixed <- sigma[fixed_vars, fixed_vars, drop=FALSE]
                sigma_free_fixed <- sigma[free_vars, fixed_vars, drop=FALSE]
                
                result$conditional_mean_offset_multiplier <-
                    sigma_free_fixed %*% solve(sigma_fixed)
            } else {
                precision_free <- precision[free_vars, free_vars,
                    drop=FALSE]
                precision_free_fixed <- precision[free_vars, fixed_vars,
                    drop=FALSE]
                
                result$conditional_mean_offset_multiplier <-
                    -1 * solve(precision_free) %*% precision_free_fixed
            }
        } else {
            result$conditional_mean_offset_multiplier <-
                conditional_mean_offset_multiplier
            if(validate_level > 1) {
                warning("did not check that provided conditional_mean_offset_multiplier agrees with sigma/precision")
            }
        }
    }
    
    return(result)
}

#' Compute mean and covariance matrix of the conditional distribution of a
#' subset of variables in a multivariate normal random vector given the other
#' variables in the vector.
#' 
#' @param x_fixed matrix of values of variables to condition on.  Each row is an
#'   observation of a subset of the variables
#' @param mean mean of the full random vector
#' @param sigma covariance matrix for the full random vector
#' @param precision precision matrix for the full random vector
#' @param fixed_vars vector of names or column/row indices corresponding to
#'   variables that we are conditioning on
#' @param free_vars vector of names or column/row indices corresponding to
#'   variables whose conditional distribution we want
#' @param conditional_sigma (OPTIONAL) conditional covariance matrix of the
#'   free variables given the fixed variables
#' @param conditional_mean_offset_multiplier (OPTIONAL) matrix used in computing
#'   the conditional mean of the free variables given the fixed variables
#' 
#' @return a named list with two entries: conditional_sigma and
#'   conditional_mean_offset_multiplier
get_conditional_mvn_params <- function(x_fixed,
    mean,
    sigma,
    precision,
    fixed_vars,
    free_vars,
    sigma_fixed,
    conditional_sigma,
    conditional_mean_offset_multiplier,
    validate_level = 1L) {
    temp <- get_conditional_mvn_intermediate_params(sigma,
        precision,
        fixed_vars,
        free_vars,
        conditional_sigma,
        conditional_mean_offset_multiplier,
        validate_level = validate_level)
    
    result <- list()
    result$conditional_sigma <- temp$conditional_sigma
    result$conditional_mean <- get_conditional_mvn_mean_from_intermediate_params(
        x_fixed = x_fixed,
        mean = mean,
        conditional_mean_offset_multiplier = temp$conditional_mean_offset_multiplier,
        fixed_vars = fixed_vars,
        free_vars = free_vars)
    
    return(result)
}

#' Compute means of the conditional distribution of a subset of variables in a
#' multivariate normal random vector given the other variables in the vector.
#' 
#' @param x_fixed matrix of values of variables to condition on.  Each row is an
#'   observation of a subset of the variables
#' @param mean mean of the full random vector
#' @param conditional_mean_offset_multiplier matrix used in computing the
#'   conditional mean of the free variables given the fixed variables
#' @param fixed_vars vector of names or column/row indices corresponding to
#'   variables that we are conditioning on
#' @param free_vars vector of names or column/row indices corresponding to
#'   variables whose conditional distribution we want
#' 
#' @return a named list with two entries: conditional_sigma and
#'   conditional_mean_offset_multiplier
get_conditional_mvn_mean_from_intermediate_params <- function(x_fixed,
    mean,
    conditional_mean_offset_multiplier,
    fixed_vars,
    free_vars) {
    if(length(fixed_vars) > 0) {
        cond_means <- matrix(rep(mean[free_vars], nrow(x_fixed)),
            nrow=nrow(x_fixed),
            byrow=TRUE)
        
        ## compute x_fixed_i - mu_fixed, stacked in a matrix
        x_fixed_differences <- sweep(x_fixed,
            2,
            mean[fixed_vars],
            `-`)
        
        ## adjust -- this step gives us a matrix where row i is the transpose of
        ## mu_d + Sigma_{free_fixed} Sigma_{fixed}^{-1} (x_fixed_i - mu_fixed)
        cond_means <- cond_means +
            x_fixed_differences %*% t(conditional_mean_offset_multiplier)
    } else {
        cond_means <- matrix(mean[free_vars],
            nrow=1,
            byrow=TRUE)
    }
    
    return(cond_means)
}


#' Compute the (log of the) constant term in the truncated multivariate normal
#' pdf that accounts for truncation:  the integral of the mvn pdf between the
#' lower and upper limits
#' 
#' @param mean mean of the truncated multivariate normal distribution
#' @param sigma covariance matrix for the truncated multivaraite normal
#'     distribution
#' @param precision precision (inverse covariance) matrix for the truncated
#'     multivariate normal distribution
#' @param lower lower truncation endpoints
#' @param upper upper truncation endpoints
#' 
#' @return log of the the integral of the mvn pdf between the
#'     lower and upper limits
compute_trunc_const_pdtmvn <- function(mean, sigma, precision, lower, upper) {
	if(missing(sigma)) {
		sigma <- solve(precision)
	}
	
	## Ideally, we would like pmvnorm to support a log=TRUE argument to give us
	## more numerical precision in high dimensional cases.  Here we compute the
	## probability and then take its logarithm.
	exp_trunc_const <- as.numeric(
		mvtnorm::pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
	)
	
	return(log(exp_trunc_const))
}

#' Compute whether each element of x is equal to an integer, up to a specified
#' tolerance level.
#' 
#' @param x numeric
#' @param tolerance numeric tolerance for comparison of integer values
#' 
#' @return logical vector of same length as x; entry i is TRUE if
#'     x[i] is within tol of as.integer(x[i])
equals_integer <- function(x, tolerance = .Machine$double.eps ^ 0.5) {
    return(sapply(x, function(x_i) {
        return(isTRUE(all.equal(x_i, as.integer(x_i))))
    }))
}

#' Compute floor(x) - 1
#' Used as default "a" function
#' 
#' @param x numeric
#' 
#' @return floor(x) - 1
floor_x_minus_1 <- function(x) {
	return(floor(x) - 1)
}

#' Calculate the Schur complement of the block of X specified by inds
#' There might be a more numerically stable way to do this.
#' 
#' @param X A square matrix
#' @param inds Indices specifying the rows/columns of the block within X
#'   to compute the Schur complement of.
#'   
#' @return The Schur complement of the submatrix of X specified by inds:
#'   A - B D^{-1} C, where A is the submatrix of X specified by inds,
#'   B is the submatrix with rows given by inds and columns not in inds,
#'   C is the submatrix with columns given by inds and rows not in inds,
#'   and D is the submatrix with rows and columns not in inds.
calc_Schur_complement <- function(X, inds) {
    if(length(inds) == ncol(X)) {
        return(X)
    } else {
        inds_complement <- which(!(seq_len(ncol(X)) %in% inds))
        X_1 <- X[inds, inds, drop=FALSE]
        X_12 <- X[inds, inds_complement, drop=FALSE]
        X_2 <- X[inds_complement, inds_complement, drop=FALSE]
        
        return(X_1 - X_12 %*% solve(X_2) %*% t(X_12))
    }
}
