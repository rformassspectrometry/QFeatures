data(se_na2)
x <- assay(se_na2)
randna  <-  rowData(se_na2)$randna

test_that("all imputation methods", {
    m <- imputeMethods()
    m <- m[m != "mixed"]
    m <- m[m != "none"]
    m <- m[m != "nbavg"] ## see next test
    for (.m in m) {
        xx <- impute_matrix(x, method = .m)
        expect_false(any(is.na(xx)))
    }
    expect_error(impute_matrix(x, method = "mixed",
                        randna = randna,
                        mnar = "min"),
                 regexp = "mar")
    expect_error(impute_matrix(x, method = "mixed",
                               randna = randna,
                               mar = "knn"),
                 regexp = "mnar")
    expect_error(impute_matrix(x, method = "mixed",
                               mnar = "min",
                               mar = "knn"),
                 regexp = "randna")
    expect_error(impute_matrix(x, method = "mixed",
                               randna = TRUE,
                               mnar = "min",
                               mar = "knn"),
                 regexp = "randna")
    mx <- impute_matrix(x, method = "mixed",
                        randna = randna,
                        mnar = "min",
                        mar = "knn")    
    expect_false(any(is.na(mx)))
})

test_that("none method", {
    xx <- impute_matrix(x, method = "none")
    expect_identical(x, xx)
})

## test_that("nbavg methods", {
##     x <- matrix(1:25, 5)
##     ## default min value
##     x[1, 2] <- 0.1
##     ## imputes as min value (or use-defined k)
##     x[1, 1] <- x[5, 5] <- NA
##     x[2, 1:2] <- NA ## [2, 1] will be min
##                     ## [2, 2] will be avg 6.05
##     ## remaing NA
##     x[3, 3:4] <- NA
##     ## average imputation
##     x[5, 2] <- NA ## will be 10
##     x[4, 3] <- NA ## will be 14
##     rownames(x) <- colnames(x) <-
##             LETTERS[1:5]

##     xx <- impute_matrix(x, "nbavg")
##     expect_true(exprs(xx[1, 2]) == 0.1)
##     expect_true(exprs(xx[1, 1]) == 0.1)
##     expect_true(exprs(xx[2, 1]) == 0.1)
##     expect_true(exprs(xx[2, 2]) == 6.05)
##     expect_true(all(is.na(exprs(xx[3, 3:4]))))
##     expect_true(exprs(xx[5, 2]) == 10)
##     expect_true(exprs(xx[4, 3]) == 14)

##     xx <- impute(x, "nbavg", k = 0)
##     expect_true(exprs(xx[1, 2]) == 0.1)
##     expect_true(exprs(xx[1, 1]) == 0)
##     expect_true(exprs(xx[2, 1]) == 0)
##     expect_true(exprs(xx[2, 2]) == 6)
##     expect_true(all(is.na(exprs(xx[3, 3:4]))))
##     expect_true(exprs(xx[5, 2]) == 10)
##     expect_true(exprs(xx[4, 3]) == 14)
## })


test_that("seed is not set by knn imputation method", {
  rand <- sapply(1:10, function(idx){
      xx <- suppressWarnings(impute_matrix(x, "knn"))
      rnorm(1)
  })
  expect_gt(max(rand) - min(rand), 0)
})
