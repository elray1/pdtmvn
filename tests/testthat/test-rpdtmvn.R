## all functions in rpdtmvn.R
##   leading X means test written,
##   leading I means implicitly tested by another test
##   leading S means simple enough that no test is required
##   no leading character means test needs to be written still
## 
##  rpdtmvn

context("rpdtmvn -- no tests written")

x_minus_0.5 <- function(x) {
    x - 0.5
}

x_plus_0.5 <- function(x) {
    x + 0.5
}

## Used as discretization function --
## consistently rounds x so that _.5 -> _+1
round_up_.5 <- function(x) {
    return(sapply(x, function(x_i) {
        if(x_i - floor(x_i) >= 0.5) {
            return(ceiling(x_i))
        } else {
            return(floor(x_i))
        }
    }))
}


#test_that("rpdtmvn_sample_w_fixed works -- 2 continuous, 1 discrete, conditioning on 1 continous and 1 discrete, truncation, sigma, upper integration limit = point + 0.5", {
#    n <- 1000
#    
#    full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
#    simulated_values <- rpdtmvn(n,
#        mean = c(va = 3, vb = 4, vc = 5),
#        sigma = full_sigma,
#        lower = c(va = -1, vb = -0.5, vc = -1.5),
#        upper = c(va = 2, vb = 1.5, vc = 3),
#        continuous_vars = 1:2,
#        discrete_vars = 3,
#        discrete_var_range_functions = list(
#            va = list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer, discretizer = round_up_.5),
#            vb = list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer, discretizer = round_up_.5),
#            vc = list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer, discretizer = round_up_.5)),
#        validate_level = 1)
#    
#    inds <- which(x[, 3] == floor(x[, 3]) & 
#            x[, 1] >= -1 & x[, 1] <= 2 &
#            x[, 2] >= -0.5 & x[, 2] <= 1.5 &
#            x[, 3] >= -1.5 & x[, 3] <= 3)
#    cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
#    cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
#    b <- floor(x[inds, 3]) + 0.5
#    a <- b - 1
#    
#    w <- matrix(NA, nrow = n, ncol = 3)
#    x_fixed <- matrix(rtmvnorm(2))
#    fixed_vars <- 2:3
#    
#    w_actual <- rpdtmvn_sample_w_fixed(n,
#        w,
#        x_fixed,
#        fixed_vars,
#        mean,
#        sigma,
#        precision,
#        lower,
#        upper,
#        sigma_continuous,
#        conditional_sigma_discrete,
#        conditional_mean_discrete_offset_multiplier,
#        continuous_vars,
#        discrete_vars,
#        discrete_var_range_functions,
#        validate_level = TRUE)
#    
#    w_expected <- ;
#    
#    expect_equal(dpdtmvn_result_log, dmvnorm_result_log, tolerance = 10^(-5))
#    expect_equal(dpdtmvn_result, dmvnorm_result, tolerance = 10^(-5))
#})


##test_that("rpdtmvn works works -- 2 continuous, 1 discrete, no conditioning, truncation, sigma, upper integration limit = point + 0.5", {
#    n <- 1000
#    
#    full_sigma <- matrix(c(4^2, 2^2, 1.5^2, 2^2, 3.5^2, 0.7^2, 1.5^2, 0.7^2, 5^2), nrow = 3, ncol = 3)
#    simulated_values <- rpdtmvn(n,
#        mean = c(va = 3, vb = 4, vc = 5),
#        sigma = full_sigma,
#        lower = c(va = -1, vb = -0.5, vc = -1.5),
#        upper = c(va = 2, vb = 1.5, vc = 3),
#        continuous_vars = 1:2,
#        discrete_vars = 3,
#        discrete_var_range_functions = list(
#            va = list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer, discretizer = round_up_.5),
#            vb = list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer, discretizer = round_up_.5),
#            vc = list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer, discretizer = round_up_.5)),
#        validate_level = 1)
#    
#    inds <- which(x[, 3] == floor(x[, 3]) & 
#            x[, 1] >= -1 & x[, 1] <= 2 &
#            x[, 2] >= -0.5 & x[, 2] <= 1.5 &
#            x[, 3] >= -1.5 & x[, 3] <= 3)
#    cond_means <- 5 + (sweep(x[inds, 1:2], 2, c(3, 4), `-`)) %*% t(full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]))
#    cond_var <- full_sigma[3, 3, drop=FALSE] - full_sigma[3, 1:2, drop=FALSE] %*% solve(full_sigma[1:2, 1:2]) %*% full_sigma[1:2, 3, drop=FALSE]
#    b <- floor(x[inds, 3]) + 0.5
#    a <- b - 1
#    
#    dmvnorm_result_log <- rep(-Inf, n)
#    dmvnorm_result_log <- rep(-Inf, n)
#    dmvnorm_result_log[inds] <- dmvnorm(x[inds, 1:2, drop=FALSE],
#            mean = 3:4,
#            sigma = full_sigma[1:2, 1:2],
#            log = TRUE) +
#        logspace_sub(pnorm(b, mean = cond_means, sd = sqrt(cond_var), log = TRUE),
#            pnorm(a, mean = cond_means, sd = sqrt(cond_var), log = TRUE)) -
#        log(pmvnorm(lower = c(-1, -0.5, -1.5),
#                upper = c(2, 1.5, 3),
#                mean = 3:5,
#                sigma = full_sigma))
#    
#    
#    dpdtmvn_result <- dpdtmvn(x,
#        mean = 3:5,
#        sigma = full_sigma,
#        lower = c(-1, -0.5, -1.5),
#        upper = c(2, 1.5, 3),
#        continuous_vars = 1:2,
#        discrete_vars = 3,
#        discrete_var_range_functions = list(list(a = x_minus_0.5, b = x_plus_0.5, in_range = equals_integer)),
#        log = FALSE,
#        validate_level = 1)
#    
#    dmvnorm_result <- exp(dmvnorm_result_log)
#    
#    expect_equal(dpdtmvn_result_log, dmvnorm_result_log, tolerance = 10^(-5))
#    expect_equal(dpdtmvn_result, dmvnorm_result, tolerance = 10^(-5))
#
#})
