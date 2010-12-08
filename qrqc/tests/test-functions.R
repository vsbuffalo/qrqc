## test-functions.R - unit test some functions

test_that(".trimRightCols", {
  t1 <- matrix(c(2, 3, 4, 5, 0, 4, 0, 0), 2, 4)
  expect_that(.trimRightCols(t1), is_equivalent_to(t1[, -4]))
})

test_that(".trimArray", {
  t1 <- c(2, 3, 4, 5, 0, 4, 0, 0)
  expect_that(.trimArray(t1), is_equivalent_to(t1[1:6]))
})


test_that(".length2weights", {
  l <- c(0, 0, 0, 0, 10, 20, 30, 40)
  expect_that(lengths2weights(l),
              is_equivalent_to(c(100, 100, 100, 100, 100, 90, 70, 40)))
})

## test_that(".meansFromBins", {
##   m <- 
##   expect_that(lengths2weights(l),
##               is_equivalent_to(c(100, 100, 100, 100, 100, 90, 70, 40)))
## })


