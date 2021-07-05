context("Test psuedo MLE clayton with test data")

test_est_cop_par <- function(data) {
  return(est_cop_par(copula = "clayton", data = data, niter = 50))
}


test_that("Test psuedo MLE on test_data has the correct result", {
  expect_equal(test_est_cop_par(test_data), 
               data.frame(outer_theta=1.376211, inner_theta=1.963197), 
               tolerance=0.1)
})
#> Test passed 🎊
