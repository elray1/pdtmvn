#' Evaluate the density of a pdTMVN distribution
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
#' @param log Logical; if TRUE, return value is log(density).
#' @param validate_level Numeric; if 0, no validation is performed (risky!).
#'   If 1, parameters are validated but warnings about checks not performed are
#'   not issued.  If 2, parameters are validated and warnings about checks not
#'   performed are issued.
#' 
#' @return Density of pdTMVN distribution evaluated at x
dpdtmvn <- function(x,
		mean,
		sigma,
		precision,
		lower = rep(-Inf, length = length(mean)),
		upper = rep(Inf, length = length(mean)),
		norm_const,
		trunc_const,
		sigma_continuous,
		conditional_sigma_discrete,
		conditional_mean_discrete_offset_multiplier,
		continuous_vars,
		discrete_vars,
		discrete_var_range_fns = lapply(seq_along(discrete_vars), function(dv) {
			list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")
		}),
		log = FALSE,
        validate_in_support = TRUE,
		validate_level = 1L) {
	
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
	
	## The rest of the code here is closely based on that in the
	## tmvtnorm package, but we
	##  - do the integration for partial discretization
	##  - allow for precomputation of trunc_const
	##    (and possibly norm_const in the future)
	
	## in_support are row indices for observations in support
    if(validate_in_support) {
        in_support <- in_pdtmvn_support(x,
            lower = lower,
            upper = upper,
            continuous_vars,
            discrete_vars,
            discrete_var_range_fns)
    } else {
        in_support <- seq_len(nrow(x))
    }
    
	## compute result -- computations are on log scale
	log_result <- rep(-Inf, nrow(x))
	
	## Compute contribution from observations of continuous variables
	if(length(continuous_vars) > 0 && length(in_support) > 0) {
		log_result[in_support] <- mvtnorm::dmvnorm(x[in_support, continuous_vars, drop=FALSE],
			mean = mean[continuous_vars],
			sigma = sigma_continuous,
			log = TRUE)
	}
	
	## Compute contribution from observations of discrete variables
	if(length(discrete_vars) > 0 && length(in_support) > 0) {
		## get conditional means of discrete vars given continuous vars
        cond_means <- get_conditional_mvn_mean_from_intermediate_params(
            x_fixed = x[in_support, continuous_vars, drop=FALSE],
            mean = mean,
            conditional_mean_offset_multiplier = conditional_mean_discrete_offset_multiplier,
            fixed_vars = continuous_vars,
            free_vars = discrete_vars)
		
		## get lower and upper bounds on integration for each discrete
        ## observation
		a_x_discrete <- plyr::laply(seq_along(discrete_vars), function(discrete_var_ind) {
			do.call(discrete_var_range_fns[[discrete_var_ind]][["a"]],
				list(x=x[in_support, discrete_vars[discrete_var_ind]])
			)
		})
        b_x_discrete <- plyr::laply(seq_along(discrete_vars), function(discrete_var_ind) {
            do.call(discrete_var_range_fns[[discrete_var_ind]][["b"]],
                list(x=x[in_support, discrete_vars[discrete_var_ind]])
            )
        })
        if(identical(length(discrete_vars), 1L)) {
            a_x_discrete <- matrix(a_x_discrete)
            b_x_discrete <- matrix(b_x_discrete)
        } else {
            a_x_discrete <- t(a_x_discrete)
            b_x_discrete <- t(b_x_discrete)
        }
        
		if(identical(length(discrete_vars), 1L)) {
			p_lt_b <- pnorm(b_x_discrete[, 1],
				mean = as.vector(cond_means),
				sd = sqrt(as.vector(conditional_sigma_discrete)),
				log = TRUE)
			
			p_lt_a <- pnorm(a_x_discrete[, 1],
				mean = as.vector(cond_means),
				sd = sqrt(as.vector(conditional_sigma_discrete)),
				log = TRUE)
			
			if(length(continuous_vars) > 0) {
				log_result[in_support] <- log_result[in_support] +
					logspace_sub(p_lt_b, p_lt_a)
			} else {
				log_result[in_support] <- logspace_sub(p_lt_b, p_lt_a)
			}
		} else {
			stop("dpdtmvn does not currently support discrete_vars with length > 1")
			## pmvnorm does not support a log=TRUE option; I suspect that we need that
			## for the probabilities to be non-zero.
			if(length(continuous_vars) > 0) {
				log_result[in_support] <- log_result[in_support] +
					apply(matrix(seq_along(in_support)), 1, function(support_row_ind) {
						mvtnorm::pmvnorm(lower=a_x_discrete[support_row_ind, ],
										upper=b_x_discrete[support_row_ind, ],
										mean=mean[discrete_vars],
										sigma=conditional_sigma_discrete)
					})
			} else {
				log_result[in_support] <- 
					apply(matrix(seq_along(in_support)), 1, function(support_row_ind) {
						mvtnorm::pmvnorm(lower=a_x_discrete[support_row_ind, ],
										upper=b_x_discrete[support_row_ind, ],
										mean=mean[discrete_vars],
										sigma=conditional_sigma_discrete)
					})
			}
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
