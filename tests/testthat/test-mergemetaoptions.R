context("test-mergemetaoptions")

test_that("merge works", {

  x <- list(a = 'a', b = NA)
  y <- list(a = NA, b = 'b')
  z <- list(a = 'a', b = 'b')
  expect_equal(mergeMetaoptions(x, y), z)

  x <- list(a = 'a', b = NA, c = NA)
  y <- list(a = NA, b = 'b')
  z <- list(a = 'a', b = 'b', c = NA)
  expect_equal(mergeMetaoptions(x, y), z)

  x <- list(a = 'a', b = NA)
  y <- list(a = NA, b = 'b', c = NA)
  z <- list(a = 'a', b = 'b', c = NA)
  expect_equal(mergeMetaoptions(x, y), z)

  x <- list(a = 'a', b = NA, c = NA)
  y <- list(a = NA, b = 'b', c = NA)
  z <- list(a = 'a', b = 'b', c = NA)
  expect_equal(mergeMetaoptions(x, y), z)

  x <- list(a = 'a', b = NA, c = NA)
  y <- list(a = 'c', b = 'b', c = NA)
  z <- list(a = 'a', b = 'b', c = NA)
  expect_warning(z2 <- mergeMetaoptions(x, y))
  expect_equal(z2, z)

  x <- list(a = 'c', b = NA, c = NA)
  y <- list(a = 'a', b = 'b', c = NA)
  z <- list(a = 'c', b = 'b', c = NA)
  expect_warning(z2 <- mergeMetaoptions(x, y))
  expect_equal(z2, z)

})
