library(pdtmvn)
library(mvtnorm)
library(tmvtnorm)

test_that("dpdtmvn works -- 1 continuous, 0 discrete, no truncation, sigma", {
	n <- 100
	x <- matrix(rnorm(n))
	
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = 5,
		sigma = matrix(4^2),
		lower = -Inf,
		upper = Inf,
		continuous_vars = 1,
		log = TRUE,
		validate_level = 1)
	
	dmvnorm_result_log <- dmvnorm(x,
		mean = 5,
		sigma = matrix(4^2),
		log = TRUE)

	dpdtmvn_result <- dpdtmvn(x,
		mean = 5,
		sigma = matrix(4^2),
		lower = -Inf,
		upper = Inf,
		continuous_vars = 1,
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- dmvnorm(x,
		mean = 5,
		sigma = matrix(4^2),
		log = FALSE)

	expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
	expect_equal(dpdtmvn_result, dmvnorm_result)
})

test_that("dpdtmvn works -- 1 continuous, 0 discrete, no truncation, precision", {
	n <- 100
	x <- matrix(rnorm(n))
	
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = 5,
		precision = matrix(1 / 4^2),
		lower = -Inf,
		upper = Inf,
		continuous_vars = 1,
		log = TRUE,
		validate_level = 1)
	
	dmvnorm_result_log <- dmvnorm(x,
		mean = 5,
		sigma = matrix(4^2),
		log = TRUE)

	dpdtmvn_result <- dpdtmvn(x,
		mean = 5,
		precision = matrix(1 / 4^2),
		lower = -Inf,
		upper = Inf,
		continuous_vars = 1,
		log = FALSE,
		validate_level = 1)
	
	dmvnorm_result <- dmvnorm(x,
		mean = 5,
		sigma = matrix(4^2),
		log = FALSE)

	expect_equal(dpdtmvn_result_log, dmvnorm_result_log)
	expect_equal(dpdtmvn_result, dmvnorm_result)
})

test_that("dpdtmvn works -- 0 continuous, 1 discrete, no truncation, sigma, some values out of discrete var range", {
	n <- 100
	x <- matrix(rnorm(n))
	x[c(2, 10:95), 1] <- floor(x[c(2, 10:95), 1])
	
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = 5,
		sigma = matrix(4^2),
		lower = -Inf,
		upper = Inf,
		discrete_vars = 1,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor")),
		log = TRUE,
		validate_level = 1)
	
	dtmvnorm_result_log <- rep(-Inf, n)
	inds <- which(as.vector(x) == floor(as.vector(x)))
	b <- floor(as.vector(x)[inds])
	a <- b - 1
	dtmvnorm_result_log[inds] <- log(pnorm(b, mean = 5, sd = 4, log = FALSE) -
		pnorm(a, mean = 5, sd = 4, log = FALSE))

	dpdtmvn_result <- dpdtmvn(x,
		mean = 5,
		sigma = matrix(4^2),
		lower = -Inf,
		upper = Inf,
		discrete_vars = 1,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor")),
		log = FALSE,
		validate_level = 1)
	
	dtmvnorm_result <- exp(dtmvnorm_result_log)

	expect_equal(dpdtmvn_result_log, dtmvnorm_result_log)
	expect_equal(dpdtmvn_result, dtmvnorm_result)
})

test_that("dpdtmvn works -- 0 continuous, 1 discrete, no truncation, precision, some values out of discrete var range", {
	n <- 100
	x <- matrix(rnorm(n))
	x[c(2, 10:95), 1] <- floor(x[c(2, 10:95), 1])
	
	dpdtmvn_result_log <- dpdtmvn(x,
		mean = 5,
		precision = matrix(1 / 4^2),
		lower = -Inf,
		upper = Inf,
		discrete_vars = 1,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor")),
		log = TRUE,
		validate_level = 1)
	
	dtmvnorm_result_log <- rep(-Inf, n)
	inds <- which(as.vector(x) == floor(as.vector(x)))
	b <- floor(as.vector(x)[inds])
	a <- b - 1
	dtmvnorm_result_log[inds] <- log(pnorm(b, mean = 5, sd = 4, log = FALSE) -
		pnorm(a, mean = 5, sd = 4, log = FALSE))

	dpdtmvn_result <- dpdtmvn(x,
		mean = 5,
		precision = matrix(1 / 4^2),
		lower = -Inf,
		upper = Inf,
		discrete_vars = 1,
		discrete_var_range_functions = list(list(a = "floor_x_minus_1", b = "floor")),
		log = FALSE,
		validate_level = 1)
	
	dtmvnorm_result <- exp(dtmvnorm_result_log)

	expect_equal(dpdtmvn_result_log, dtmvnorm_result_log)
	expect_equal(dpdtmvn_result, dtmvnorm_result)
})

