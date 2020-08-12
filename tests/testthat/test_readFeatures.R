data(hlpsms)
f <- tempfile()
x <- hlpsms
write.csv(x, file = f, row.names = FALSE)


test_that("readQFeatures", {
    ft1 <- readQFeatures(x, ecol = 1:10, name = "psms")
    ft2 <- readQFeatures(f, ecol = 1:10, name = "psms")
    expect_equal(ft1, ft2)
    ft1 <- readQFeatures(x, ecol = 1:10, name = NULL)
    ft2 <- readQFeatures(f, ecol = 1:10, name = NULL)
    expect_equal(ft1, ft2)
    expect_message(ft1 <- readQFeatures(x, ecol = 1:10, name = "psms", fname = "Sequence"),
                   "Making assay rownames unique.")
    ft2 <- readQFeatures(f, ecol = 1:10, name = "psms", fname = "Sequence")
    expect_equal(ft1, ft2)
    ft2 <- readQFeatures(x, ecol = 1:10, name = "psms", fname = 11)
    ft3 <- readQFeatures(f, ecol = 1:10, name = "psms", fname = 11)
    expect_equal(ft1, ft2)
    expect_equal(ft1, ft3)
    ecol <- c("X126", "X127C", "X127N", "X128C", "X128N", "X129C",
              "X129N", "X130C", "X130N", "X131")
    ft1 <- readQFeatures(x, ecol = ecol, name = "psms")
    ft2 <- readQFeatures(f, ecol = ecol, name = "psms")
    expect_equal(ft1, ft2)
    ecol <- LETTERS[1:10]
    expect_error(readQFeatures(x, ecol = ecol, name = "psms"))
    expect_error(readQFeatures(f, ecol = ecol, name = "psms"))
    expect_error(readQFeatures(x, ecol = 1:10, name = "psms",
                              fname = "not_present"))
    expect_error(readQFeatures(f, ecol = 1:10, name = "psms",
                              fname = "not_present"))
})

test_that("readQFeatures", {
    ft1 <- readSummarizedExperiment(x, ecol = 1:10)
    ft2 <- readSummarizedExperiment(f, ecol = 1:10)
    expect_equal(ft1, ft2)
    ft1 <- readSummarizedExperiment(x, ecol = 1:10, fname = "Sequence")
    ft2 <- readSummarizedExperiment(f, ecol = 1:10, fname = "Sequence")
    expect_equal(ft1, ft2)
    ft3 <- readSummarizedExperiment(f, ecol = 1:10, fname = 11)
    expect_equal(ft1, ft3)
    ecol <- LETTERS[1:10]
    expect_error(readSummarizedExperiment(x, ecol = ecol, name = "psms"))
    expect_error(readSummarizedExperiment(f, ecol = ecol, name = "psms"))
    expect_error(readSummarizedExperiment(x, ecol = 1:10, name = "psms",
                                          fname = "not_present"))
    expect_error(readSummarizedExperiment(f, ecol = 1:10, name = "psms",
                                          fname = "not_present"))
    expect_true(inherits(ft1, "SummarizedExperiment"))
})
