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

test_that("dpdtmvn works -- 0 continuous, 1 discrete", {
	;
})

test_that("dpdtmvn works -- 2 continuous, 1 discrete", {
	;
})

