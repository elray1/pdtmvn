### numerical functions

## interface to R's C API for logspace arithmetic

logspace_sub <- function(logx, logy) {
	if(!identical(length(logx), length(logy))) {
		stop("logx and logy must have the same length")
	} else {
#		return(sapply(seq_along(logx), function(ind) {
#			return(.Call("logspace_sub_C", as.numeric(logx[ind]), as.numeric(logy[ind])))
#		}))
        return(logspace_sub_matrix_rows(cbind(logx, logy)))
	}
}

logspace_add <- function(logx, logy) {
	if(!identical(length(logx), length(logy))) {
		stop("logx and logy must have the same length")
	} else {
		return(sapply(seq_along(logx), function(ind) {
			return(.Call("logspace_add_C", as.numeric(logx), as.numeric(logy)))
		}))
	}
}

logspace_sum_matrix_rows <- function(logX) {
	return(.Call("logspace_sum_matrix_rows_C", as.numeric(logX), as.integer(nrow(logX)), as.integer(ncol(logX))))
}

logspace_sub_matrix_rows <- function(logX) {
	if(!is.matrix(logX) || !identical(ncol(logX), 2L))
		stop("logX must be a matrix with 2 columns")

	return(.Call("logspace_sub_matrix_rows_C", as.numeric(logX), as.integer(nrow(logX))))
}
