library(pdtmvn)
library(microbenchmark)

lower <- c(-Inf, 0.7, -10)
upper <- c(Inf, 1, 5)

x <- rbind(
	c(-1, 5, 3), # upper bound exceeded on second entry
	c(0, 0, 1.3), # lower bound exceeded on second entry
	c(-Inf, 1, 0), # infinite
	c(Inf, 1, 0), # infinite
	c(0, 1, 5) # passes
)

in_trunc_supp_1 <- function(x, lower, upper) {
	apply(x, 1, function(x_row) {
		all(x_row >= lower & x_row <= upper & !any(is.infinite(x_row)))
	})
}

in_trunc_supp_2 <- function(x, lower, upper) {
	rowSums(
		sweep(x, 2, lower, "<") +
		sweep(x, 2, upper, ">") +
		is.infinite(x)
	) == 0
}

in_trunc_supp_3 <- function(x, lower, upper) {
	apply(
		sweep(x, 2, lower, ">=") &
		sweep(x, 2, upper, "<=") &
		is.finite(x),
		1,
		all
	)
}

in_trunc_supp_4 <- function(x, lower, upper) {
	apply(x, 1, function(x_row) {
		all(x_row >= lower & x_row <= upper & is.finite(x_row))
	})
}

in_trunc_supp_5 <- function(x, lower, upper) {
	rowSums(
		sweep(x, 2, lower, "<") |
		sweep(x, 2, upper, ">") |
		is.infinite(x)
	) == 0
}

in_trunc_supp_6 <- function(x, lower, upper) {
	apply(x, 1, function(x_row) {
		for(i in seq_along(x_row)) {
			if(x_row[i] < lower || x_row[i] > upper || is.infinite(x_row[i])) {
				return(FALSE)
			}
		}
		return(TRUE)
	})
}

in_trunc_supp_7 <- function(x, lower, upper) {
	apply(x, 1, function(x_row) {
		any(x_row < lower | x_row > upper | is.infinite(x_row))
	})
}

in_trunc_supp_8 <- function(x, lower, upper) {
	storage.mode(x) <- "double"
	storage.mode(lower) <- "double"
	storage.mode(upper) <- "double"
	
	return(as.logical(.Call("in_truncation_domain_C", x, lower, upper)))
}


for(fn_num in c(2:5, 8)) {
	fn_name <- paste0("in_trunc_supp_", fn_num)
	print(identical(
		in_trunc_supp_1(x, lower, upper),
		do.call(fn_name, list(x = x, lower = lower, upper = upper))
	))
	stopifnot(identical(
		in_trunc_supp_1(x, lower, upper),
		do.call(fn_name, list(x = x, lower = lower, upper = upper))
	))
}

microbenchmark(
	in_trunc_supp_1(x, lower, upper),
	in_trunc_supp_4(x, lower, upper),
	in_trunc_supp_7(x, lower, upper),
	in_trunc_supp_2(x, lower, upper),
	in_trunc_supp_3(x, lower, upper),
	in_trunc_supp_5(x, lower, upper),
	in_trunc_supp_6(x, lower, upper),
	in_trunc_supp_8(x, lower, upper),
	times = 10000L
)


