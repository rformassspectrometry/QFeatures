data(hlpsms)
x <- hlpsms

test_that("readQFeatures", {
    ft1 <- readQFeatures(x, ecol = 1:10, name = "psms")
    ft2 <- readQFeatures(x, colAnnotation = 1:10, name = "psms")
    expect_equal(ft1, ft2)
    ft1 <- readQFeatures(x, ecol = 1:10, name = NULL)
    ft2 <- readQFeatures(x, colAnnotation = 1:10, name = NULL)
    expect_equal(ft1, ft2)
    expect_message(ft1 <- readQFeatures(x, ecol = 1:10, name = "psms", fnames = "Sequence"),
                   "Making assay rownames unique.")
    ft2 <- readQFeatures(x, colAnnotation = 1:10, name = "psms", fnames = "Sequence")
    expect_equal(ft1, ft2)
    ft2 <- readQFeatures(x, ecol = 1:10, name = "psms", fname = 11)
    ft3 <- readQFeatures(x, colAnnotation = 1:10, name = "psms", fname = 11)
    expect_equal(ft1, ft2)
    expect_equal(ft1, ft3)
    ecol <- c("X126", "X127C", "X127N", "X128C", "X128N", "X129C",
              "X129N", "X130C", "X130N", "X131")
    ft1 <- readQFeatures(x, ecol = ecol, name = "psms")
    ft2 <- readQFeatures(x, colAnnotation = ecol, name = "psms")
    expect_equal(ft1, ft2)
    ecol <- LETTERS[1:10]
    expect_error(readQFeatures(x, ecol = ecol, name = "psms"))
    expect_error(readQFeatures(x, colAnnotation = ecol, name = "psms"))
    expect_error(readQFeatures(x, ecol = 1:10, name = "psms",
                               fnames = "not_present"))
    expect_error(readQFeatures(x, ecol = 1:10, name = "psms",
                               fnames = "not_present"))
})

test_that("readSummarizedExperiment", {
    ft1 <- readSummarizedExperiment(x, ecol = 1:10)
    ft2 <- readSummarizedExperiment(x, colAnnotation = 1:10)
    expect_equal(ft1, ft2)
    ft3 <- readSummarizedExperiment(x, ecol = 1:10, fnames = "Sequence")
    ft4 <- readSummarizedExperiment(x, colAnnotation = 1:10, fnames = "Sequence")
    expect_equal(ft3, ft4)
    ft5 <- readSummarizedExperiment(x, ecol = 1:10, fname = 11)
    expect_equal(ft3, ft5)
    ## Read data with only 1 quantitation column
    ft5 <- readSummarizedExperiment(x, ecol = 1, fname = 11)
    ## Check column names
    ecol <- c("X126", "X127C", "X127N", "X128C", "X128N", "X129C",
              "X129N", "X130C", "X130N", "X131")
    expect_identical(colnames(ft3), ecol)
    expect_identical(colnames(ft5), ecol[1])
    ## Provide ecol as logical
    ecol <- seq_along(x) %in% 1:10
    expect_identical(ft1, readSummarizedExperiment(x, ecol = ecol))

    ## Expect errors
    ecol <- LETTERS[1:10]
    expect_error(readSummarizedExperiment(x, ecol = ecol, name = "psms"))
    expect_error(readSummarizedExperiment(x, colAnnotation = ecol, name = "psms"))
    expect_error(readSummarizedExperiment(x, ecol = 1:10, name = "psms",
                                          fnames = "not_present"))
    expect_error(readSummarizedExperiment(f, colAnnotation = 1:10, name = "psms",
                                          fnames = "not_present"))
    expect_true(inherits(ft1, "SummarizedExperiment"))
})
