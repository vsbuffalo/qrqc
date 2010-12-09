## test-functions.R - unit test some functions
EPSILON <- 0.001
  
test_that(".trimRightCols", {
  t1 <- matrix(c(2, 3, 4, 5, 0, 4, 0, 0), 2, 4)
  expect_that(.trimRightCols(t1), is_equivalent_to(t1[, -4]))
})

test_that(".trimArray", {
  t1 <- c(2, 3, 4, 5, 0, 4, 0, 0)
  expect_that(.trimArray(t1), is_equivalent_to(t1[1:6]))
})


test_that("length2weights", {
  l <- c(0, 0, 0, 0, 10, 20, 30, 40)
  expect_that(lengths2weights(l),
              is_equivalent_to(c(100, 100, 100, 100, 100, 90, 70, 40)))
})

test_that("meanFromBins", {
  m1 <- structure(list(`26` = c(2L, 0L, 0L, 2L, 3L, 0L, 0L, 0L, 1L, 0L),
                       `27` = c(1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L),
                       `28` = c(3L, 7L, 5L, 4L, 4L, 1L, 0L, 3L, 2L, 2L),
                       `29` = c(1L, 0L, 1L, 1L, 1L, 3L, 0L, 4L, 1L, 1L),
                       `30` = c(4L, 2L, 5L, 3L, 4L, 3L, 3L, 4L, 3L, 5L),
                       `31` = c(0L, 0L, 2L, 0L, 1L, 1L, 1L, 1L, 1L, 3L), 
                       `32` = c(3L, 8L, 3L, 2L, 7L, 5L, 0L, 4L, 8L, 4L),
                       `33` = c(4L, 1L, 3L, 3L, 4L, 2L, 1L, 5L, 1L, 5L),
                       `34` = c(2L, 6L, 3L, 1L, 3L, 2L, 12L, 4L, 4L, 1L),
                       `35` = c(10L, 6L, 11L, 14L, 8L, 13L, 8L, 13L, 12L, 13L),
                       `36` = c(14L, 13L, 16L, 15L, 14L, 18L, 21L, 15L, 20L, 16L)),
                  .Names = c("26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"),
                  row.names = c(NA, 10L), class = "data.frame")
  m1 <- cbind(position=1:nrow(m1), m1)
  expect_that(meanFromBins(m1) - 33.60487 < EPSILON, is_true())
})

test_that("binned2quantilefunc", {
  b <- c(0, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 0)
  names(b) <- 1:length(b)
  f <- binned2quantilefunc(b)
  expect_that(f(0.25), is_equivalent_to(5.45))
  expect_that(f(0.5), is_equivalent_to(7.5))
  expect_that(f(0.75), is_equivalent_to(9.55))
})

test_that("binned2boxplot", {
  b <- c(0, 1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 0)
  names(b) <- 1:length(b)
  f <- binned2quantilefunc(b)
  ans <- structure(c(2, 3.63333333333333, 5.45, 7.5, 9.55, 11.3666666666667, 14),
                   .Names = c("ymin", "alt.lower", "lower", "middle", "upper", "alt.upper", "ymax"))
                                                  
  expect_that(binned2boxplot(b), is_equivalent_to(ans))
})
