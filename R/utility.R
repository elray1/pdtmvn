### numerical functions

## interface to R's C API for logspace arithmetic


#' An interface to R's C function logspace_add
#' Computes log(exp(log_x) + exp(log_y))
#' Unlike R's implementation, we return
#' -INFINITY if both log_x and log_y are -INFINITY
#' 
#' @param log_x a numeric vector
#' @param log_y a numeric vector with the same length as log_x
#' 
#' @return a numeric vector with the same length as log_x and log_y with the
#'   computed value
logspace_sub <- function(logx, logy) {
	if(!identical(length(logx), length(logy))) {
		stop("logx and logy must have the same length")
	} else {
        return(logspace_sub_matrix_rows(cbind(logx, logy)))
	}
}

#' An interface to R's C function logspace_add
#' Computes log(exp(log_x) + exp(log_y))
#' Unlike R's implementation, we return
#' -INFINITY if both log_x and log_y are -INFINITY
#' This could be sped up by calling logspace_sum_matrix_rows with a similar
#' pattern as in logspace_sub
#' 
#' @param logx a numeric vector
#' @param logy a numeric vector with the same length as logx
#' 
#' @return a numeric vector with the same length as logx and logy with the
#'   computed value
logspace_add <- function(logx, logy) {
	if(!identical(length(logx), length(logy))) {
		stop("logx and logy must have the same length")
	} else {
		return(sapply(seq_along(logx), function(ind) {
			return(.Call("logspace_add_C", as.numeric(logx), as.numeric(logy)))
		}))
	}
}

#' Given a numeric matrix logX, compute the row-sums of exp(logX) in log space
#' 
#' @param logX is a numeric matrix object
#' 
#' @return a numeric vector of length nrow(logX) with the row sums of exp(logX)
logspace_sum_matrix_rows <- function(logX) {
	return(.Call("logspace_sum_matrix_rows_C", as.numeric(logX), as.integer(nrow(logX)), as.integer(ncol(logX))))
}

#' Given a numeric matrix logX with two columns, compute
#'   exp(logX[, 1]) - exp(logX[, 2]) in log space
#' 
#' @param logX is a numeric matrix object with two columns
#' 
#' @return a numeric vector of length nrow(logX) with the row sums of exp(logX)
logspace_sub_matrix_rows <- function(logX) {
	if(!is.matrix(logX) || !identical(ncol(logX), 2L))
		stop("logX must be a matrix with 2 columns")

	return(.Call("logspace_sub_matrix_rows_C", as.numeric(logX), as.integer(nrow(logX))))
}
