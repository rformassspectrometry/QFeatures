data("ft_na")
se_na <- ft_na[["na"]]
## Create data matrices containing 0's or NA's
minf <- m0 <- m <- assay(se_na)
m0[3, 2] <- m0[4, 1] <- m0[1, 1] <- 0L
minf[m0 == 0] <- Inf
## Create SE objects containing 0's or NA's
se_zero <- se_inf <- se_na
assay(se_zero) <- m0
assay(se_inf) <- minf
## Create a QFeatures object containing 0's or NA's
ft0 <- QFeatures(list(na = se_na, zero = se_zero, inf = se_inf),
                 colData = DataFrame(row.names = LETTERS[1:3]))

test_that("function: .zeroIsNA, .infIsNA", {
    ## .zeroIsNA
    expect_identical(se_na, QFeatures:::.zeroIsNA(se_zero))
    expect_identical(se_na, zeroIsNA(se_zero))
    ## .infIsNA
    ## Note: use expect_equal because m is int and minf is num
    expect_equal(se_na, QFeatures:::.infIsNA(se_inf))
    expect_equal(se_na, infIsNA(se_inf))
})

test_that("function: .nNAByAssay, .nNAByMargin, .nNA, and .nNAi", {
    ## .nNAByAssay
    nNAassay <- QFeatures:::.nNAByAssay(se_na)
    ## The expected results are initialized after manual inspection
    expect_identical(nNAassay,
                     DataFrame(nNA = 3L,
                               pNA = 3 / 12))

    ## .nNAByMargin
    nNArows <- QFeatures:::.nNAByMargin(se_na, MARGIN = 1)
    nNAcols <- QFeatures:::.nNAByMargin(se_na, MARGIN = 2)
    ## The expected results are initialized after manual inspection
    expect_identical(nNArows,
                     DataFrame(name = rownames(se_na),
                               nNA = c(1L, 0L, 1L, 1L),
                               pNA = c(1/3, 0, 1/3, 1/3)))
    expect_identical(nNAcols,
                     DataFrame(name = colnames(se_na),
                               nNA = c(2L, 1L, 0L),
                               pNA = c(1/2, 1/4, 0)))

    ## .nNA for SummarizedExperiemnt
    expect_identical(QFeatures:::.nNA(se_na),
                     list(nNA = nNAassay, nNArows = nNArows,
                          nNAcols =nNAcols))
    ## Expect only 0's (no missing data) for se_zero
    expect_true(all(sapply(QFeatures:::.nNA(se_zero),
                           function(x) all(x[, "pNA"] == 0))))

    ## .nNAi for QFeatures
    ## The expected results are initialized after manual inspection
    nNAassay <- c(3L, 0L)
    pNAassay <- nNAassay / 12
    nNArows <- c(1L, 0L, 1L, 1L, rep(0L, 4))
    pNArows <- nNArows / 3
    nNAcols <- c(2L, 1L, 0L, rep(0L, 3))
    pNAcols <- nNAcols / 4
    ## Test results
    n_na <- QFeatures:::.nNAi(ft0, 1:2)
    ## .nNAByAssay
    expect_identical(n_na$nNA,
                     DataFrame(assay = names(ft0)[1:2],
                               nNA = nNAassay, pNA = pNAassay))
    ## .nNAByMargin by row
    expect_identical(n_na$nNArows,
                     DataFrame(assay = rep(names(ft0)[1:2], each = 4),
                               name = unlist(rownames(ft0)[1:2],
                                             use.names = FALSE),
                               nNA = nNArows, pNA = pNArows))
    ## .nNAByMargin by column
    expect_identical(n_na$nNAcols,
                     DataFrame(assay = rep(names(ft0)[1:2], each = 3),
                               name = unlist(colnames(ft0)[1:2],
                                             use.names = FALSE),
                               nNA = nNAcols,
                               pNA = pNAcols))
    ## Check .nNAi with character indexing
    expect_identical(n_na, QFeatures:::.nNAi(ft0, c("na", "zero")))
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
    ## pNA must be in [0, 1]
    expect_error(QFeatures:::.row_for_filterNA(assay(se_na), pNA = -1))
    expect_error(QFeatures:::.row_for_filterNA(assay(se_na), pNA = 2))
})

test_that("zeroIsNA,QFeatures", {
    ## Check that no "zeroIsNA,missing" exist
    expect_error(ft <- zeroIsNA(ft0))
    ## Check replace on 0 containing assay
    ft <- zeroIsNA(ft0, "zero")
    expect_identical(se_na, ft[["zero"]])
    ## Check replace on NA containing assay (no replacement)
    ft <- zeroIsNA(ft0, "na")
    expect_identical(se_na, ft[["na"]])
    ## zeroIsNA on multiple assays
    ft <- zeroIsNA(ft0, c(1L, 2L)) ## test zeroIsNA,integer
    expect_identical(ft, zeroIsNA(ft, c("na", "zero"))) ## test zeroIsNA,character
    expect_identical(ft, zeroIsNA(ft, c(1.1, 2.1))) ## test zeroIsNA,numeric
    expect_identical(assay(ft[["na"]]), assay(se_na))
    expect_identical(assay(ft[["zero"]]), assay(se_na))
})

test_that("infIsNA,QFeatures", {
    ## Check that no "infIsNA,missing" exist
    expect_error(ft <- infIsNA(ft0))
    ## Check replace on inf containing assay
    ft <- infIsNA(ft0, "inf")
    ## Note: use expect_equal because m is int and minf is num
    expect_equal(se_na, ft[["inf"]])
    ## Check replace on NA containing assay (no replacement)
    ft <- infIsNA(ft0, "na")
    expect_equal(se_na, ft[["na"]])
    ## infIsNA on multiple assays
    ft <- infIsNA(ft0, c(1L, 3L)) ## test infIsNA,integer
    expect_identical(ft, infIsNA(ft, c("na", "inf"))) ## test infIsNA,character
    expect_identical(ft, infIsNA(ft, c(1, 3))) ## test infIsNA,numeric
    expect_equal(assay(ft[["inf"]]), assay(se_na))
    expect_equal(assay(ft[["inf"]]), assay(se_na))
})

test_that("nNA,SummarizedExperiment and nNA,QFeatures", {
    ## Add an assay with different dimensions (cf issue 118)
    ft0 <- addAssay(ft0, ft0[[1]][1:2, 1:2], name = "subset1")
    ## Method vs internal function
    expect_identical(nNA(se_na), QFeatures:::.nNA(se_na))
    expect_identical(nNA(ft0, 1:4), QFeatures:::.nNAi(ft0, 1:4))
    ## nNA on a single assay
    expect_identical(nNA(ft0[[1]]), nNA(se_na))
    ## nNA on multiple assays
    expect_identical(nNA(ft0, 1:3), nNA(ft0, c("na", "zero", "inf")))
})

test_that("filterNA,QFeatures and filterNA,SummarizedExperiment", {
    ## filterNA on SE
    ## se_na contains 1 NA in 3 rows. Removing all rows with NA results
    ## in an assay with 1 row and 3 colums
    se_na_filtered <- filterNA(se_na, pNA = 0)
    expect_identical(dim(se_na_filtered), c(1L, 3L))
    expect_identical(assay(se_na_filtered), m[2, , drop = FALSE])
    ## filterNA on QFeatures
    ft_filtered <- filterNA(ft0, i = seq_along(ft0), pNA = 0)
    expect_identical(se_na_filtered, ft_filtered[[1]])
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
    ## Check filterNA on linked assays
    expect_true(
        validObject(filterNA(ft_na, i = seq_along(ft_na), pNA = 0.5))
    )
})
