test_that("topIdx throws errors", {
    m <- matrix(1:(6*3), ncol = 3, byrow = TRUE)
    i <- rep(1:3, each = 2)
    expect_error(
        QFeatures:::topIdx(m, i, n = 0L),
        "must be an integer >= 1.")
    expect_error(
        QFeatures:::topIdx(m, i, n = -1L),
        "must be an integer >= 1.")
    expect_error(
        QFeatures:::topIdx(m, i, n = NA),
        "must be an integer >= 1.")
    expect_error(
        QFeatures:::topIdx(m, 1:4, n = 2L),
        "must be equal.")
})

test_that("topIdx works", {
    m <- matrix(1:(6*3), ncol = 3, byrow = TRUE)
    ## groups are ((1, 2)(3, 4)(5, 6))
    i <- rep(1:3, each = 2)
    ## second of each group
    expect_identical(
        c(2L, 4L, 6L),
        QFeatures:::topIdx(m, i, n = 1L,
                           fun = rowSums,
                           decreasing = TRUE))
    ## first of each group
    expect_identical(
        c(1L, 3L, 5L),
        QFeatures:::topIdx(m, i, n = 1L,
                           fun = rowSums,
                           decreasing = FALSE))

    ## groups are ((1, 2, 3)(4, 5, 6))
    i <- rep(1:2, each = 3)
    ## last and second of each group
    expect_identical(
        c(3L, 2L, 6L, 5L),
        QFeatures:::topIdx(m, i, n = 2L,
                           fun = rowSums,
                           decreasing = TRUE))
    expect_identical(
        c(1L, 2L, 4L, 5L),
        QFeatures:::topIdx(m, i, n = 2L,
                           fun = rowSums,
                           decreasing = FALSE))
    ## all from each group, from smaller to largest
    expect_identical(
        1:6,
        QFeatures:::topIdx(m, i, n = 3L,
                           fun = rowSums,
                           decreasing = FALSE))
    ## all from each group, from largest to smallest
    expect_identical(
        c(3:1,6:4),
        QFeatures:::topIdx(m, i, n = 3L,
                           fun = rowSums,
                           decreasing = TRUE))

})


test_that("topIdx works with NA", {
    m <- matrix(1:(6*3), ncol = 3, byrow = TRUE)
    m[4, 3] <- NA
    m[4, 2] <- 100
    ## groups are ((1, 2)(3, 4)(5, 6))
    i <- rep(1:3, each = 2)
    ## row 4 drops due to NA
    expect_identical(
        c(2L, 3L, 6L),
        QFeatures:::topIdx(m, i, n = 1L,
                           fun = rowSums,
                           decreasing = TRUE))
    ## recover row 4 when ignoring NA
    expect_identical(
        c(2L, 4L, 6L),
        QFeatures:::topIdx(m, i, n = 1L,
                           fun = rowSums,
                           decreasing = TRUE,
                           na.rm = TRUE))
})

test_that("filterTopFeatures,SummarizedExperiment works", {
    data(feat1)
    se <- feat1[[1]]
    res <- filterTopFeatures(se, n = 2L, fcol = "Sequence")
    ## order of peptides (in Sequence) are ELGNDAYK (PSM4, 5, 6),
    ## IAEESNFPFIK (PSM7, 8, 9, 10) and SYGFNAAR (PSM1, 2, 3). Among
    ## these the last and penultimate are chose, in that order.
    expect_identical(res, se[c(6, 5, 10, 9, 3, 2), ])
})

test_that("filterTopFeatures,SummarizedExperiment throws errors", {
    data(feat1)
    se <- feat1[[1]]
    expect_error(
        res <- filterTopFeatures(se, n = 2L),
        "'fcol' not found in the assay's rowData.")
    expect_error(
        res <- filterTopFeatures(se, n = 2L, fcol = "---"),
        "'fcol' not found in the assay's rowData.")
})
