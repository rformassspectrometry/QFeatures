data(ft_na)
se_na <- ft_na[["na"]]
minf <- m0 <- m <- assay(se_na)
m0[3, 2] <- m0[4, 1] <- m0[1, 1] <- 0
minf[m0 == 0] <- Inf
se_zero <- SummarizedExperiment(assay = m0)
se_inf <- SummarizedExperiment(assay = minf)

ft0 <- QFeatures(list(na = se_na, zero = se_zero, inf = se_inf),
                 colData = DataFrame(row.names = LETTERS[1:3]))


test_that("function: .zeroIsNA, .infIsNA, .nNA, and .nNAi", {
    ## .zeroIsNA
    expect_equivalent(se_na, QFeatures:::.zeroIsNA(se_zero))
    expect_equivalent(se_na, zeroIsNA(se_zero))
    ## .infIsNA
    expect_equivalent(se_na, QFeatures:::.infIsNA(se_inf))
    expect_equivalent(se_na, infIsNA(se_inf))
    ## .nNA
    n_na <- QFeatures:::.nNA(se_na)
    expect_identical(n_na[[1]], 3/(3 * 4))
    expect_identical(n_na[[2]], table(c(0, 1, 1, 1)))
    expect_identical(n_na[[3]], c(A = 2, B = 1, C = 0))
    expect_identical(nNA(se_na), QFeatures:::.nNA(se_na))
    n_na <- QFeatures:::.nNA(se_zero)
    expect_identical(n_na[[1]], 0)
    expect_identical(n_na[[2]], table(c(0, 0, 0, 0)))
    expect_identical(n_na[[3]], c(A = 0, B = 0, C = 0))
    expect_identical(nNA(se_zero), QFeatures:::.nNA(se_zero))
    ## .nNAi
    n_na <- QFeatures:::.nNAi(ft0, 1:2)
    expect_identical(n_na, QFeatures:::.nNAi(ft0, c("na", "zero")))
    expect_identical(n_na, QFeatures:::.nNAi(ft0, c(1.1, 2.1)))
    expect_identical(n_na[[1]], c(na = 0.25, zero = 0.00))
    expect_identical(n_na[[2]],
                     matrix(c(1, 4, 3, rep(0, 5)), nrow = 2,
                            dimnames = list(c("na", "zero"), 0:3)))
    expect_identical(n_na[[3]],
                     matrix(c(2, 1, rep(0, 4)), nrow = 2, byrow = TRUE,
                            dimnames = list(c("na", "zero"), LETTERS[1:3])))
    expect_identical(nNA(ft0, 1:2), QFeatures:::.nNAi(ft0, 1:2))
})


test_that("function: .row_for_filterNA", {
    def <- QFeatures:::.row_for_filterNA(m)
    def_0 <- QFeatures:::.row_for_filterNA(m, pNA = 0L)
    expect_error(QFeatures:::.row_for_filterNA(se_na))
    expect_error(QFeatures:::.row_for_filterNA(m, pNA = TRUE))
    expect_error(QFeatures:::.row_for_filterNA(m, pNA = "0"))
    expect_error(QFeatures:::.row_for_filterNA(m, pNA = c(A = 0, B = 0.5, C = 1)))
    expect_identical(def, def_0)
    expect_identical(def, c(a = FALSE, b = TRUE, c = FALSE, d = FALSE))
    expect_identical(QFeatures:::.row_for_filterNA(assay(se_zero)),
                     c(a = TRUE, b = TRUE, c = TRUE, d = TRUE))
    expect_identical(QFeatures:::.row_for_filterNA(assay(se_na), pNA = .9),
                     c(a = TRUE, b = TRUE, c = TRUE, d = TRUE))
    expect_identical(QFeatures:::.row_for_filterNA(assay(se_zero)),
                     c(a = TRUE, b = TRUE, c = TRUE, d = TRUE))
    expect_identical(QFeatures:::.row_for_filterNA(assay(se_na), pNA = .5),
                     c(a = TRUE, b = TRUE, c = TRUE, d = TRUE))
    expect_identical(QFeatures:::.row_for_filterNA(assay(se_na), pNA = .33),
                     QFeatures:::.row_for_filterNA(assay(se_na), pNA = 0))
    expect_identical(QFeatures:::.row_for_filterNA(assay(se_na), pNA = 0),
                     QFeatures:::.row_for_filterNA(assay(se_na), pNA = -1))
    expect_identical(QFeatures:::.row_for_filterNA(assay(se_na), pNA = 1),
                     QFeatures:::.row_for_filterNA(assay(se_na), pNA = 2))
})


test_that("zeroIsNA,QFeatures", {
    expect_error(ft <- zeroIsNA(ft0))
    ft <- zeroIsNA(ft0, 1)
    expect_equivalent(ft[["na"]], ft[["zero"]])
    ft <- zeroIsNA(ft0, "na")
    expect_equivalent(ft[["na"]], ft[["zero"]])
    ## zeroIsNA on multiple assays
    ft <- zeroIsNA(ft0, 1:2)
    expect_identical(ft, zeroIsNA(ft, c("na", "zero")))
    expect_identical(ft, zeroIsNA(ft, c(1.1, 2.1)))
    expect_equivalent(assay(ft[["na"]]), assay(se_na))
    expect_equivalent(assay(ft[["zero"]]), assay(se_na))
})

test_that("infIsNA,QFeatures", {
    expect_error(ft <- infIsNA(ft0))
    ft <- zeroIsNA(ft0, 1)
    expect_equivalent(ft[["na"]], ft[["zero"]])
    ft <- zeroIsNA(ft0, "na")
    expect_equivalent(ft[["na"]], ft[["zero"]])
    ## zeroIsNA on multiple assays
    ft <- zeroIsNA(ft0, 1:2)
    expect_identical(ft, zeroIsNA(ft, c("na", "zero")))
    expect_identical(ft, zeroIsNA(ft, c(1.1, 2.1)))
    expect_equivalent(assay(ft[["na"]]), assay(se_na))
    expect_equivalent(assay(ft[["zero"]]), assay(se_na))
})


test_that("nNA,QFeatures", {
    expect_error(nNA(ft0))
    n_na <- nNA(ft0, i = 1:2)
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
    expect_identical(nNA(ft0[[1]]), nNA(ft0, i = 1))
    ## nNA on multiple assays
    n_na <- nNA(ft0, 1:2)
    expect_identical(n_na, nNA(ft0, c("na", "zero")))
    expect_identical(n_na[[1]], c(na = 0.25, zero = 0.00))
    expect_identical(n_na[[2]],
                     matrix(c(1, 4, 3, rep(0, 5)), nrow = 2,
                            dimnames = list(c("na", "zero"), 0:3)))
    expect_identical(n_na[[3]],
                     matrix(c(2, 1, rep(0, 4)), nrow = 2, byrow = TRUE,
                            dimnames = list(c("na", "zero"), LETTERS[1:3])))
})


test_that("filterNA,QFeatures and filterNA,SummarizedExperiment", {
    se_na_filtered <- filterNA(se_na)
    expect_error(filterNA(ft0))
    ft_filtered <- filterNA(ft0, i = seq_along(ft0))
    expect_equivalent(se_na_filtered, ft_filtered[[1]])
    expect_identical(assay(se_na_filtered), m[2, , drop = FALSE])
    se_na_filtered <- filterNA(se_na, pNA = 0.9)
    ft_filtered <- filterNA(ft0, i = seq_along(ft0), pNA = 0.9)
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
                    Y = LETTERS[1:2],
                    .n = c(2L, 2L),
                    row.names = 1:2)
    expect_equivalent(rowData(ft_na[[2]]), rd)
})
