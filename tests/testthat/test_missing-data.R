data(ft_na)
se_na <- ft_na[["na"]]
m0 <- m <- assay(se_na)
m0[3, 2] <- m0[4, 1] <- m0[1, 1] <- 0
se_zero <- SummarizedExperiment(assay = m0)
               
ft0 <- Features(list(na = se_na, zero = se_zero),
                colData = DataFrame(row.names = LETTERS[1:3]))


test_that("function: .zeroIsNA and .nNA", {
    expect_equivalent(se_na, Features:::.zeroIsNA(se_zero))
    expect_equivalent(se_na, zeroIsNA(se_zero))    
    n_na <- Features:::.nNA(se_na)
    expect_identical(n_na[[1]], 3/(3 * 4))
    expect_identical(n_na[[2]], table(c(0, 1, 1, 1)))
    expect_identical(n_na[[3]], c(A = 2, B = 1, C = 0))    
    expect_identical(nNA(se_na), Features:::.nNA(se_na))
    n_na <- Features:::.nNA(se_zero)    
    expect_identical(n_na[[1]], 0)
    expect_identical(n_na[[2]], table(c(0, 0, 0, 0)))
    expect_identical(n_na[[3]], c(A = 0, B = 0, C = 0))
    expect_identical(nNA(se_zero), Features:::.nNA(se_zero))
})



test_that("function: .row_for_filterNA", {
    def <- Features:::.row_for_filterNA(m)
    def_0 <- Features:::.row_for_filterNA(m, pNA = 0L)
    expect_error(Features:::.row_for_filterNA(se_na))
    expect_error(Features:::.row_for_filterNA(m, pNA = TRUE))
    expect_error(Features:::.row_for_filterNA(m, pNA = "0"))
    expect_error(Features:::.row_for_filterNA(m, pNA = c(A = 0, B = 0.5, C = 1)))
    expect_identical(def, def_0)
    expect_identical(def, c(a = FALSE, b = TRUE, c = FALSE, d = FALSE)) 
    expect_identical(Features:::.row_for_filterNA(assay(se_zero)),
                     c(a = TRUE, b = TRUE, c = TRUE, d = TRUE))
    expect_identical(Features:::.row_for_filterNA(assay(se_na), pNA = .9),
                     c(a = TRUE, b = TRUE, c = TRUE, d = TRUE))                    
    expect_identical(Features:::.row_for_filterNA(assay(se_zero)),
                     c(a = TRUE, b = TRUE, c = TRUE, d = TRUE))
    expect_identical(Features:::.row_for_filterNA(assay(se_na), pNA = .5),
                     c(a = TRUE, b = TRUE, c = TRUE, d = TRUE))                     
    expect_identical(Features:::.row_for_filterNA(assay(se_na), pNA = .33),
                     Features:::.row_for_filterNA(assay(se_na), pNA = 0))
    expect_identical(Features:::.row_for_filterNA(assay(se_na), pNA = 0),
                     Features:::.row_for_filterNA(assay(se_na), pNA = -1))
    expect_identical(Features:::.row_for_filterNA(assay(se_na), pNA = 1),
                     Features:::.row_for_filterNA(assay(se_na), pNA = 2))
})

test_that("function: zeroIsNA and nNA", {
    n_na <- Features:::.nNA(se_na)
    expect_identical(n_na[[1]], 3/(3 * 4))
    expect_identical(n_na[[2]], table(c(c(0, 1, 1, 1))))
    expect_identical(n_na[[3]], c(A = 2, B = 1, C = 0))
})


test_that("zeroIsNA,Features", {
    ft <- zeroIsNA(ft0)
    expect_equivalent(ft[["na"]], ft[["zero"]])
    ft <- zeroIsNA(ft0, 1)
    expect_equivalent(ft[["na"]], ft[["zero"]])
    ft <- zeroIsNA(ft0, "na")
    expect_equivalent(ft[["na"]], ft[["zero"]])
})

test_that("nNA,Features", {
    n_na <- nNA(ft0)
    expect_identical(n_na[[1]], c(na = 3/(4*3), zero = 0))
    expect_identical(n_na[[2]],
                     matrix(c(1, 3, 0, 0, 4, 0, 0, 0),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(c("na", "zero"),
                                            0:3)))
    expect_identical(n_na[[3]],
                     matrix(c(2, 1, 0, 0, 0, 0),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(c("na", "zero"),
                                            LETTERS[1:3])))
    expect_identical(nNA(ft0, 1L), nNA(se_na))
    expect_identical(nNA(ft0, 1), nNA(se_na))
    expect_identical(nNA(ft0, "na"), nNA(se_na))
    expect_identical(nNA(ft0[, , 1]), nNA(ft0[[1]]))
})


test_that("filterNA,Features and filterNA,SummarizedExperiment", {
    se_na_filtered <- filterNA(se_na)
    ft_filtered <- filterNA(ft0)
    expect_equivalent(se_na_filtered, ft_filtered[[1]])
    expect_identical(assay(se_na_filtered), m[2, , drop = FALSE])
    se_na_filtered <- filterNA(se_na, pNA = 0.9)
    ft_filtered <- filterNA(ft0, pNA = 0.9)
    expect_equivalent(se_na_filtered, ft_filtered[[1]])
    expect_equivalent(se_na_filtered, ft0[[2]])
})

test_that("aggregateFeatures with missing data", {
    expect_message(ft_na <- aggregateFeatures(ft_na, "na", fcol = "X", name = "agg_na",
                                              fun = colSums))
    expect_message(ft_na <- aggregateFeatures(ft_na, "na", fcol = "X", name = "agg_na_rm",
                                              fun = colSums,
                                              na.rm = TRUE))
    agg1 <- matrix(c(NA, NA, 20,
                     NA, 14, 22),
                   ncol = 3, byrow = TRUE,
                   dimnames = list(1:2, LETTERS[1:3]))
    agg2 <- matrix(c(3, 5, 20,
                     2, 14, 22),
                   ncol = 3, byrow = TRUE,
                   dimnames = list(1:2, LETTERS[1:3]))
    expect_identical(assay(ft_na[[2]]), agg1)
    expect_identical(assay(ft_na[[3]]), agg2)
    expect_identical(rowData(ft_na[[2]]), rowData(ft_na[[3]]))
    rd <- DataFrame(X = c(1L, 2L),
                    .n = c(2L, 2L),
                    row.names = 1:2)
    expect_equivalent(rowData(ft_na[[2]]), rd)
})
