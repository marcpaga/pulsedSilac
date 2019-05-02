context("test-rbindlinkerdf")

test_that("rbindlinkerdf works, no overlap", {

  linkM1 <- buildLinkerDf(protIDs = c(1:3),
                          pepIDs = c(1:4),
                          protToPep = list('1' = c(1), '2' = c(2),
                                           '3' = c(3, 4)))

  linkM2 <- buildLinkerDf(protIDs = c(4:6),
                          pepIDs = c(5:8),
                          protToPep = list('4' = c(5), '5' = c(6),
                                           '6' = c(7, 8)))

  new.lm <- buildLinkerDf(protIDs = c(1:6),
                          pepIDs = c(1:8),
                          protToPep = c(list('1' = c(1), '2' = c(2),
                                             '3' = c(3, 4)),
                                        list('4' = c(5), '5' = c(6),
                                             '6' = c(7, 8))))

  expect_silent(lm <- rbind_linkerDf(linkM1, linkM2))
  expect_equal(lm, new.lm)
})

test_that("rbindlinkerdf works, protein overlap", {

  linkM1 <- buildLinkerDf(protIDs = c(1:3),
                          pepIDs = c(1:4),
                          protToPep = list('1' = c(1), '2' = c(2),
                                           '3' = c(3, 4)))

  linkM2 <- buildLinkerDf(protIDs = c(3:5),
                          pepIDs = c(5:8),
                          protToPep = list('3' = c(5), '4' = c(6),
                                           '5' = c(7, 8)))

  new.lm <- buildLinkerDf(protIDs = c(1:5),
                          pepIDs = c(1:8),
                          protToPep = c(list('1' = c(1), '2' = c(2),
                                             '3' = c(3, 4, 5)),
                                        list('4' = c(6), '5' = c(7, 8))))

  expect_silent(lm <- rbind_linkerDf(linkM1, linkM2))
  expect_equal(lm, new.lm)
})


test_that("rbindlinkerdf works, peptide overlap", {

  linkM1 <- buildLinkerDf(protIDs = c(1:3),
                          pepIDs = c(1:4),
                          protToPep = list('1' = c(1), '2' = c(2),
                                           '3' = c(3, 4)))

  linkM2 <- buildLinkerDf(protIDs = c(4:7),
                          pepIDs = c(3:8),
                          protToPep = list('4' = c(3), '5' = c(4),
                                           '6' = c(5, 6), '7' = c(7, 8)))

  new.lm <- buildLinkerDf(protIDs = c(1:7),
                          pepIDs = c(1:8),
                          protToPep = c(list('1' = c(1), '2' = c(2),
                                             '3' = c(3, 4)),
                                        list('4' = c(3), '5' = c(4),
                                             '6' = c(5, 6), '7' = c(7, 8))))

  expect_silent(lm <- rbind_linkerDf(linkM1, linkM2))
  expect_equal(lm, new.lm)
})

test_that("rbindlinkerdf works, both overlap", {

  linkM1 <- buildLinkerDf(protIDs = c(1:3),
                          pepIDs = c(1:4),
                          protToPep = list('1' = c(1),
                                           '2' = c(2),
                                           '3' = c(3, 4)))

  linkM2 <- buildLinkerDf(protIDs = c(3:7),
                          pepIDs = c(3:8),
                          protToPep = list('3' = c(3, 6, 8),
                                           '4' = c(3),
                                           '5' = c(4),
                                           '6' = c(5, 6),
                                           '7' = c(7, 8)))

  new.lm <- buildLinkerDf(protIDs = c(1:7),
                          pepIDs = c(1:8),
                          protToPep = c(list('1' = c(1),
                                             '2' = c(2),
                                             '3' = c(3, 4, 6, 8)),
                                        list('4' = c(3),
                                             '5' = c(4),
                                             '6' = c(5, 6),
                                             '7' = c(7, 8))))

  expect_silent(lm <- rbind_linkerDf(linkM1, linkM2))
  expect_equal(lm, new.lm)
})
