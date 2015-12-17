library(pdtmvn)
library(mvtnorm)
library(tmvtnorm)

################################################################################
## Tests with upper integration limit = point
################################################################################

test_that("dpdtmvn works -- 1 continuous, 1 discrete, no truncation, sigma, upper integration limit = point", {
	n <- 100
	x <- matrix(rnorm(2*n), nrow = n)
	x[c(2, 10:95), 2] <- floor(x[c(2, 10:95), 1])
	
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = c(3, 5),
		sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2),
		lower = rep(-Inf, 2),
		upper = rep(Inf, 2),
		continuous_vars = 1,
		discrete_vars = 2,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = TRUE,
		validate_level = 1)
	
	inds <- which(x[, 2] == floor(x[, 2]))
	cond_means <- 5 + (2^2) / (4^2) * (x[inds, 1] - 3)
	cond_var <- 5^2 - 2^2 * (1 / 4^2) * 2^2
	b <- floor(as.vector(x[inds, 2]))
	a <- b - 1
	
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1, drop=FALSE],
		mean = 3,
		sigma = matrix(4^2),
		log = TRUE) +
		log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
			pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE))
		

	dpdtmvn_result <- dpdtmvn(x,
		mean = c(3, 5),
		sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2),
		lower = rep(-Inf, 2),
		upper = rep(Inf, 2),
		continuous_vars = 1,
		discrete_vars = 2,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- exp(dmvnorm_result_log)

	expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
	expect_equal(dpdtmvn_result, dmvnorm_result)
})

test_that("dpdtmvn works -- 1 continuous, 1 discrete, no truncation, precision, upper integration limit = point", {
	n <- 100
	x <- matrix(rnorm(2*n), nrow = n)
	x[c(2, 10:95), 2] <- floor(x[c(2, 10:95), 1])
	
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = c(3, 5),
		precision = solve(matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)),
		lower = rep(-Inf, 2),
		upper = rep(Inf, 2),
		continuous_vars = 1,
		discrete_vars = 2,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = TRUE,
		validate_level = 1)
	
	inds <- which(x[, 2] == floor(x[, 2]))
	cond_means <- 5 + (2^2) / (4^2) * (x[inds, 1] - 3)
	cond_var <- 5^2 - 2^2 * (1 / 4^2) * 2^2
	b <- floor(as.vector(x[inds, 2]))
	a <- b - 1
	
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1, drop=FALSE],
		mean = 3,
		sigma = matrix(4^2),
		log = TRUE) +
		log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
		pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE))
		

	dpdtmvn_result <- dpdtmvn(x,
		mean = c(3, 5),
		precision = solve(matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)),
		lower = rep(-Inf, 2),
		upper = rep(Inf, 2),
		continuous_vars = 1,
		discrete_vars = 2,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- exp(dmvnorm_result_log)

	expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
	expect_equal(dpdtmvn_result, dmvnorm_result)
})



test_that("dpdtmvn works -- 1 continuous, 1 discrete, truncation, sigma, upper integration limit = point", {
	n <- 100
	x <- matrix(rnorm(2*n), nrow = n)
	x[c(2, 10:95), 2] <- floor(x[c(2, 10:95), 1])
	
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = c(3, 5),
		sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2),
		lower = c(-1, -0.5),
		upper = c(2, 1.5),
		continuous_vars = 1,
		discrete_vars = 2,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = TRUE,
		validate_level = 1)
	
	inds <- which(x[, 2] == floor(x[, 2]) & 
			x[, 1] >= -1 & x[, 1] <= 2 &
			x[, 2] >= -0.5 & x[, 2] <= 1.5)
	cond_means <- 5 + (2^2) / (4^2) * (x[inds, 1] - 3)
	cond_var <- 5^2 - 2^2 * (1 / 4^2) * 2^2
	b <- floor(as.vector(x[inds, 2]))
	a <- b - 1
	
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1, drop=FALSE],
		mean = 3,
		sigma = matrix(4^2),
		log = TRUE) +
		log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
		pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE)) -
		log(pmvnorm(lower = c(-1, -0.5),
			upper = c(2, 1.5),
			mean = c(3, 5),
			sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)))
	
	
	dpdtmvn_result <- dpdtmvn(x,
		mean = c(3, 5),
		sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2),
		lower = c(-1, -0.5),
		upper = c(2, 1.5),
		continuous_vars = 1,
		discrete_vars = 2,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- exp(dmvnorm_result_log)

	expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
	expect_equal(dpdtmvn_result, dmvnorm_result)
})


test_that("dpdtmvn works -- 1 continuous, 1 discrete, truncation, precision, upper integration limit = point", {
	n <- 100
	x <- matrix(rnorm(2*n), nrow = n)
	x[c(2, 10:95), 2] <- floor(x[c(2, 10:95), 1])
	
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = c(3, 5),
		precision = solve(matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)),
		lower = c(-1, -0.5),
		upper = c(2, 1.5),
		continuous_vars = 1,
		discrete_vars = 2,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = TRUE,
		validate_level = 1)
	
	inds <- which(x[, 2] == floor(x[, 2]) & 
			x[, 1] >= -1 & x[, 1] <= 2 &
			x[, 2] >= -0.5 & x[, 2] <= 1.5)
	cond_means <- 5 + (2^2) / (4^2) * (x[inds, 1] - 3)
	cond_var <- 5^2 - 2^2 * (1 / 4^2) * 2^2
	b <- floor(as.vector(x[inds, 2]))
	a <- b - 1
	
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1, drop=FALSE],
		mean = 3,
		sigma = matrix(4^2),
		log = TRUE) +
		log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
		pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE)) -
		log(pmvnorm(lower = c(-1, -0.5),
			upper = c(2, 1.5),
			mean = c(3, 5),
			sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)))
	
	
	dpdtmvn_result <- dpdtmvn(x,
		mean = c(3, 5),
		precision = solve(matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)),
		lower = c(-1, -0.5),
		upper = c(2, 1.5),
		continuous_vars = 1,
		discrete_vars = 2,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- exp(dmvnorm_result_log)

	expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
	expect_equal(dpdtmvn_result, dmvnorm_result)
})



test_that("dpdtmvn works -- 2 continuous, 1 discrete, no truncation, sigma, upper integration limit = point", {
	n <- 100
	x <- matrix(rnorm(3*n), nrow = n)
	x[c(2, 10:95), 3] <- floor(x[c(3, 10:95), 1])
	
	full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = 3:5,
		sigma = full_sigma,
		lower = rep(-Inf, 3),
		upper = rep(Inf, 3),
		continuous_vars = 1:2,
		discrete_vars = 3,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = TRUE,
		validate_level = 1)
	
	inds <- which(x[, 3] == floor(x[, 3]))
	cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
	cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
	b <- floor(x[inds, 3])
	a <- b - 1
	
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1:2, drop=FALSE],
		mean = 3:4,
		sigma = full_sigma[1:2, 1:2],
		log = TRUE) +
		log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
		pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE))
		

	dpdtmvn_result <- dpdtmvn(x,
		mean = 3:5,
		sigma = full_sigma,
		lower = rep(-Inf, 3),
		upper = rep(Inf, 3),
		continuous_vars = 1:2,
		discrete_vars = 3,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- exp(dmvnorm_result_log)

	expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
	expect_equal(dpdtmvn_result, dmvnorm_result)
})

test_that("dpdtmvn works -- 2 continuous, 1 discrete, no truncation, precision, upper integration limit = point", {
	n <- 100
	x <- matrix(rnorm(3*n), nrow = n)
	x[c(2, 10:95), 3] <- floor(x[c(3, 10:95), 1])
	
	full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = 3:5,
		precision = solve(full_sigma),
		lower = rep(-Inf, 3),
		upper = rep(Inf, 3),
		continuous_vars = 1:2,
		discrete_vars = 3,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = TRUE,
		validate_level = 1)
	
	inds <- which(x[, 3] == floor(x[, 3]))
	cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
	cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
	b <- floor(x[inds, 3])
	a <- b - 1
	
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1:2, drop=FALSE],
		mean = 3:4,
		sigma = full_sigma[1:2, 1:2],
		log = TRUE) +
		log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
		pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE))
		

	dpdtmvn_result <- dpdtmvn(x,
		mean = 3:5,
		precision = solve(full_sigma),
		lower = rep(-Inf, 3),
		upper = rep(Inf, 3),
		continuous_vars = 1:2,
		discrete_vars = 3,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- exp(dmvnorm_result_log)

	expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
	expect_equal(dpdtmvn_result, dmvnorm_result)
})



test_that("dpdtmvn works -- 2 continuous, 1 discrete, truncation, sigma, upper integration limit = point", {
	n <- 100
	x <- matrix(rnorm(3*n), nrow = n)
	x[c(2, 10:95), 3] <- floor(x[c(3, 10:95), 1])
	
	full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = 3:5,
		sigma = full_sigma,
		lower = c(-1, -0.5, -1.5),
		upper = c(2, 1.5, 3),
		continuous_vars = 1:2,
		discrete_vars = 3,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = TRUE,
		validate_level = 1)
	
	inds <- which(x[, 3] == floor(x[, 3]) & 
			x[, 1] >= -1 & x[, 1] <= 2 &
			x[, 2] >= -0.5 & x[, 2] <= 1.5 &
			x[, 3] >= -1.5 & x[, 3] <= 3)
	cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
	cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
	b <- floor(x[inds, 3])
	a <- b - 1
	
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1:2, drop=FALSE],
		mean = 3:4,
		sigma = full_sigma[1:2, 1:2],
		log = TRUE) +
		logspace_sub(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = TRUE),
			pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = TRUE)) -
		log(pmvnorm(lower = c(-1, -0.5, -1.5),
			upper = c(2, 1.5, 3),
			mean = 3:5,
			sigma = full_sigma))
	
	
	dpdtmvn_result <- dpdtmvn(x,
		mean = 3:5,
		sigma = full_sigma,
		lower = c(-1, -0.5, -1.5),
		upper = c(2, 1.5, 3),
		continuous_vars = 1:2,
		discrete_vars = 3,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- exp(dmvnorm_result_log)
	
	expect_equal(dpdtmvn_result_log, dmvnorm_result_log, tolerance = 10^(-5))
	expect_equal(dpdtmvn_result, dmvnorm_result, tolerance = 10^(-5))
})


test_that("dpdtmvn works -- 2 continuous, 1 discrete, truncation, precision, upper integration limit = point", {
	n <- 100
	x <- matrix(rnorm(3*n), nrow = n)
	x[c(2, 10:95), 3] <- floor(x[c(3, 10:95), 1])
	
	full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = 3:5,
		precision = solve(full_sigma),
		lower = c(-1, -0.5, -1.5),
		upper = c(2, 1.5, 3),
		continuous_vars = 1:2,
		discrete_vars = 3,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = TRUE,
		validate_level = 1)
	
	inds <- which(x[, 3] == floor(x[, 3]) & 
			x[, 1] >= -1 & x[, 1] <= 2 &
			x[, 2] >= -0.5 & x[, 2] <= 1.5 &
			x[, 3] >= -1.5 & x[, 3] <= 3)
	cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
	cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
	b <- floor(x[inds, 3])
	a <- b - 1
	
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log <- rep(-Inf, n)
	dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1:2, drop=FALSE],
		mean = 3:4,
		sigma = full_sigma[1:2, 1:2],
		log = TRUE) +
		logspace_sub(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = TRUE),
			pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = TRUE)) -
		log(pmvnorm(lower = c(-1, -0.5, -1.5),
			upper = c(2, 1.5, 3),
			mean = 3:5,
			sigma = full_sigma))
	
	
	dpdtmvn_result <- dpdtmvn(x,
		mean = 3:5,
		precision = solve(full_sigma),
		lower = c(-1, -0.5, -1.5),
		upper = c(2, 1.5, 3),
		continuous_vars = 1:2,
		discrete_vars = 3,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor", in_range = "equals_integer")),
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- exp(dmvnorm_result_log)
	
	expect_equal(dpdtmvn_result_log, dmvnorm_result_log, tolerance = 10^(-5))
	expect_equal(dpdtmvn_result, dmvnorm_result, tolerance = 10^(-5))
})





################################################################################
## Tests with upper integration limit = point + 0.5
################################################################################

x_minus_0.5 <- function(x) {
	x - 0.5
}

x_plus_0.5 <- function(x) {
	x + 0.5
}


test_that("dpdtmvn works -- 1 continuous, 1 discrete, no truncation, sigma, upper integration limit = point", {
    n <- 100
    x <- matrix(rnorm(2*n), nrow = n)
    x[c(2, 10:95), 2] <- floor(x[c(2, 10:95), 1])
    
    dpdtmvn_result_log <- dpdtmvn(x,
        mean = c(3, 5),
        sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2),
        lower = rep(-Inf, 2),
        upper = rep(Inf, 2),
        continuous_vars = 1,
        discrete_vars = 2,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = TRUE,
        validate_level = 1)
    
    inds <- which(x[, 2] == floor(x[, 2]))
    cond_means <- 5 + (2^2) / (4^2) * (x[inds, 1] - 3)
    cond_var <- 5^2 - 2^2 * (1 / 4^2) * 2^2
    b <- floor(as.vector(x[inds, 2])) + 0.5
    a <- b - 1
    
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1, drop=FALSE],
        mean = 3,
        sigma = matrix(4^2),
        log = TRUE) +
        log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
                pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE))
    
    
    dpdtmvn_result <- dpdtmvn(x,
        mean = c(3, 5),
        sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2),
        lower = rep(-Inf, 2),
        upper = rep(Inf, 2),
        continuous_vars = 1,
        discrete_vars = 2,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = FALSE,
        validate_level = 1)
    
    dmvnorm_result <- exp(dmvnorm_result_log)
    
    expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
    expect_equal(dpdtmvn_result, dmvnorm_result)
})

test_that("dpdtmvn works -- 1 continuous, 1 discrete, no truncation, precision, upper integration limit = point", {
    n <- 100
    x <- matrix(rnorm(2*n), nrow = n)
    x[c(2, 10:95), 2] <- floor(x[c(2, 10:95), 1])
    
    dpdtmvn_result_log <- dpdtmvn(x,
        mean = c(3, 5),
        precision = solve(matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)),
        lower = rep(-Inf, 2),
        upper = rep(Inf, 2),
        continuous_vars = 1,
        discrete_vars = 2,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = TRUE,
        validate_level = 1)
    
    inds <- which(x[, 2] == floor(x[, 2]))
    cond_means <- 5 + (2^2) / (4^2) * (x[inds, 1] - 3)
    cond_var <- 5^2 - 2^2 * (1 / 4^2) * 2^2
    b <- floor(as.vector(x[inds, 2])) + 0.5
    a <- b - 1
    
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1, drop=FALSE],
        mean = 3,
        sigma = matrix(4^2),
        log = TRUE) +
        log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
                pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE))
    
    
    dpdtmvn_result <- dpdtmvn(x,
        mean = c(3, 5),
        precision = solve(matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)),
        lower = rep(-Inf, 2),
        upper = rep(Inf, 2),
        continuous_vars = 1,
        discrete_vars = 2,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = FALSE,
        validate_level = 1)
    
    dmvnorm_result <- exp(dmvnorm_result_log)
    
    expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
    expect_equal(dpdtmvn_result, dmvnorm_result)
})



test_that("dpdtmvn works -- 1 continuous, 1 discrete, truncation, sigma, upper integration limit = point", {
    n <- 100
    x <- matrix(rnorm(2*n), nrow = n)
    x[c(2, 10:95), 2] <- floor(x[c(2, 10:95), 1])
    
    dpdtmvn_result_log <- dpdtmvn(x,
        mean = c(3, 5),
        sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2),
        lower = c(-1, -0.5),
        upper = c(2, 1.5),
        continuous_vars = 1,
        discrete_vars = 2,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = TRUE,
        validate_level = 1)
    
    inds <- which(x[, 2] == floor(x[, 2]) & 
            x[, 1] >= -1 & x[, 1] <= 2 &
            x[, 2] >= -0.5 & x[, 2] <= 1.5)
    cond_means <- 5 + (2^2) / (4^2) * (x[inds, 1] - 3)
    cond_var <- 5^2 - 2^2 * (1 / 4^2) * 2^2
    b <- floor(as.vector(x[inds, 2])) + 0.5
    a <- b - 1
    
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1, drop=FALSE],
        mean = 3,
        sigma = matrix(4^2),
        log = TRUE) +
        log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
                pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE)) -
        log(pmvnorm(lower = c(-1, -0.5),
            upper = c(2, 1.5),
            mean = c(3, 5),
            sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)))
    
    
    dpdtmvn_result <- dpdtmvn(x,
        mean = c(3, 5),
        sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2),
        lower = c(-1, -0.5),
        upper = c(2, 1.5),
        continuous_vars = 1,
        discrete_vars = 2,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = FALSE,
        validate_level = 1)
    
    dmvnorm_result <- exp(dmvnorm_result_log)
    
    expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
    expect_equal(dpdtmvn_result, dmvnorm_result)
})


test_that("dpdtmvn works -- 1 continuous, 1 discrete, truncation, precision, upper integration limit = point", {
    n <- 100
    x <- matrix(rnorm(2*n), nrow = n)
    x[c(2, 10:95), 2] <- floor(x[c(2, 10:95), 1])
    
    dpdtmvn_result_log <- dpdtmvn(x,
        mean = c(3, 5),
        precision = solve(matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)),
        lower = c(-1, -0.5),
        upper = c(2, 1.5),
        continuous_vars = 1,
        discrete_vars = 2,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = TRUE,
        validate_level = 1)
    
    inds <- which(x[, 2] == floor(x[, 2]) & 
            x[, 1] >= -1 & x[, 1] <= 2 &
            x[, 2] >= -0.5 & x[, 2] <= 1.5)
    cond_means <- 5 + (2^2) / (4^2) * (x[inds, 1] - 3)
    cond_var <- 5^2 - 2^2 * (1 / 4^2) * 2^2
    b <- floor(as.vector(x[inds, 2])) + 0.5
    a <- b - 1
    
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1, drop=FALSE],
        mean = 3,
        sigma = matrix(4^2),
        log = TRUE) +
        log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
                pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE)) -
        log(pmvnorm(lower = c(-1, -0.5),
            upper = c(2, 1.5),
            mean = c(3, 5),
            sigma = matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)))
    
    
    dpdtmvn_result <- dpdtmvn(x,
        mean = c(3, 5),
        precision = solve(matrix(c(4^2, 2^2, 2^2, 5^2), nrow = 2, ncol = 2)),
        lower = c(-1, -0.5),
        upper = c(2, 1.5),
        continuous_vars = 1,
        discrete_vars = 2,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = FALSE,
        validate_level = 1)
    
    dmvnorm_result <- exp(dmvnorm_result_log)
    
    expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
    expect_equal(dpdtmvn_result, dmvnorm_result)
})



test_that("dpdtmvn works -- 2 continuous, 1 discrete, no truncation, sigma, upper integration limit = point", {
    n <- 100
    x <- matrix(rnorm(3*n), nrow = n)
    x[c(2, 10:95), 3] <- floor(x[c(3, 10:95), 1])
    
    full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
    dpdtmvn_result_log <- dpdtmvn(x,
        mean = 3:5,
        sigma = full_sigma,
        lower = rep(-Inf, 3),
        upper = rep(Inf, 3),
        continuous_vars = 1:2,
        discrete_vars = 3,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = TRUE,
        validate_level = 1)
    
    inds <- which(x[, 3] == floor(x[, 3]))
    cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
    cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
    b <- floor(x[inds, 3]) + 0.5
    a <- b - 1
    
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1:2, drop=FALSE],
        mean = 3:4,
        sigma = full_sigma[1:2, 1:2],
        log = TRUE) +
        log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
                pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE))
    
    
    dpdtmvn_result <- dpdtmvn(x,
        mean = 3:5,
        sigma = full_sigma,
        lower = rep(-Inf, 3),
        upper = rep(Inf, 3),
        continuous_vars = 1:2,
        discrete_vars = 3,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = FALSE,
        validate_level = 1)
    
    dmvnorm_result <- exp(dmvnorm_result_log)
    
    expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
    expect_equal(dpdtmvn_result, dmvnorm_result)
})

test_that("dpdtmvn works -- 2 continuous, 1 discrete, no truncation, precision, upper integration limit = point", {
    n <- 100
    x <- matrix(rnorm(3*n), nrow = n)
    x[c(2, 10:95), 3] <- floor(x[c(3, 10:95), 1])
    
    full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
    dpdtmvn_result_log <- dpdtmvn(x,
        mean = 3:5,
        precision = solve(full_sigma),
        lower = rep(-Inf, 3),
        upper = rep(Inf, 3),
        continuous_vars = 1:2,
        discrete_vars = 3,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = TRUE,
        validate_level = 1)
    
    inds <- which(x[, 3] == floor(x[, 3]))
    cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
    cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
    b <- floor(x[inds, 3]) + 0.5
    a <- b - 1
    
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1:2, drop=FALSE],
        mean = 3:4,
        sigma = full_sigma[1:2, 1:2],
        log = TRUE) +
        log(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = FALSE) -
                pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = FALSE))
    
    
    dpdtmvn_result <- dpdtmvn(x,
        mean = 3:5,
        precision = solve(full_sigma),
        lower = rep(-Inf, 3),
        upper = rep(Inf, 3),
        continuous_vars = 1:2,
        discrete_vars = 3,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = FALSE,
        validate_level = 1)
    
    dmvnorm_result <- exp(dmvnorm_result_log)
    
    expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
    expect_equal(dpdtmvn_result, dmvnorm_result)
})



test_that("dpdtmvn works -- 2 continuous, 1 discrete, truncation, sigma, upper integration limit = point", {
    n <- 100
    x <- matrix(rnorm(3*n), nrow = n)
    x[c(2, 10:95), 3] <- floor(x[c(3, 10:95), 1])
    
    full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
    dpdtmvn_result_log <- dpdtmvn(x,
        mean = 3:5,
        sigma = full_sigma,
        lower = c(-1, -0.5, -1.5),
        upper = c(2, 1.5, 3),
        continuous_vars = 1:2,
        discrete_vars = 3,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = TRUE,
        validate_level = 1)
    
    inds <- which(x[, 3] == floor(x[, 3]) & 
            x[, 1] >= -1 & x[, 1] <= 2 &
            x[, 2] >= -0.5 & x[, 2] <= 1.5 &
            x[, 3] >= -1.5 & x[, 3] <= 3)
    cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
    cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
    b <- floor(x[inds, 3]) + 0.5
    a <- b - 1
    
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1:2, drop=FALSE],
        mean = 3:4,
        sigma = full_sigma[1:2, 1:2],
        log = TRUE) +
        logspace_sub(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = TRUE),
            pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = TRUE)) -
        log(pmvnorm(lower = c(-1, -0.5, -1.5),
            upper = c(2, 1.5, 3),
            mean = 3:5,
            sigma = full_sigma))
    
    
    dpdtmvn_result <- dpdtmvn(x,
        mean = 3:5,
        sigma = full_sigma,
        lower = c(-1, -0.5, -1.5),
        upper = c(2, 1.5, 3),
        continuous_vars = 1:2,
        discrete_vars = 3,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = FALSE,
        validate_level = 1)
    
    dmvnorm_result <- exp(dmvnorm_result_log)
    
    expect_equal(dpdtmvn_result_log, dmvnorm_result_log, tolerance = 10^(-5))
    expect_equal(dpdtmvn_result, dmvnorm_result, tolerance = 10^(-5))
})


test_that("dpdtmvn works -- 2 continuous, 1 discrete, truncation, precision, upper integration limit = point", {
    n <- 100
    x <- matrix(rnorm(3*n), nrow = n)
    x[c(2, 10:95), 3] <- floor(x[c(3, 10:95), 1])
    
    full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
    dpdtmvn_result_log <- dpdtmvn(x,
        mean = 3:5,
        precision = solve(full_sigma),
        lower = c(-1, -0.5, -1.5),
        upper = c(2, 1.5, 3),
        continuous_vars = 1:2,
        discrete_vars = 3,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = TRUE,
        validate_level = 1)
    
    inds <- which(x[, 3] == floor(x[, 3]) & 
            x[, 1] >= -1 & x[, 1] <= 2 &
            x[, 2] >= -0.5 & x[, 2] <= 1.5 &
            x[, 3] >= -1.5 & x[, 3] <= 3)
    cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
    cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
    b <- floor(x[inds, 3]) + 0.5
    a <- b - 1
    
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log <- rep(-Inf, n)
    dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1:2, drop=FALSE],
        mean = 3:4,
        sigma = full_sigma[1:2, 1:2],
        log = TRUE) +
        logspace_sub(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = TRUE),
            pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = TRUE)) -
        log(pmvnorm(lower = c(-1, -0.5, -1.5),
            upper = c(2, 1.5, 3),
            mean = 3:5,
            sigma = full_sigma))
    
    
    dpdtmvn_result <- dpdtmvn(x,
        mean = 3:5,
        precision = solve(full_sigma),
        lower = c(-1, -0.5, -1.5),
        upper = c(2, 1.5, 3),
        continuous_vars = 1:2,
        discrete_vars = 3,
        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
        log = FALSE,
        validate_level = 1)
    
    dmvnorm_result <- exp(dmvnorm_result_log)
    
    expect_equal(dpdtmvn_result_log, dmvnorm_result_log, tolerance = 10^(-5))
    expect_equal(dpdtmvn_result, dmvnorm_result, tolerance = 10^(-5))
})


