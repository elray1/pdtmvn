## all functions in pdtmvn-misc.R
##   leading X means test written,
##   leading I means implicitly tested by another test
##   leading S means simple enough that no test is required
##   no leading character means test needs to be written still
## 
## I in_pdtmvn_support
## X calc_in_truncation_support
## I validate_params_pdtmvn
## I compute_sigma_subcomponents
## I get_conditional_mvn_intermediate_params
## S get_conditional_mvn_params
## I get_conditional_mvn_mean_from_intermediate_params
## S compute_trunc_const_pdtmvn
## S equals_integer
## S floor_x_minus_1
## I calc_Schur_complement

context("pdtmvn-misc -- no tests implemented")

## We should implement tests here for numerical stability?
#test_that("calc_Schur_complement works", {
#    
#})

test_that("calc_in_truncation_support works", {
        lower <- c(-Inf, 0.7, -10)
        upper <- c(Inf, 1, 5)
        
        x <- rbind(
            c(-1, 5, 3), # upper bound exceeded on second entry
            c(0, 0, 1.3), # lower bound exceeded on second entry
            c(-Inf, 1, 0), # infinite
            c(Inf, 1, 0), # infinite
            c(0, 1, 5) # passes
        )
        
        expect_identical(
            c(FALSE, FALSE, FALSE, FALSE, TRUE),
            calc_in_truncation_support(x, lower, upper)
        )
})
