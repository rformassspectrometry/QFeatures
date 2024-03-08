data(hlpsms)
x <- hlpsms
x$file <- paste0("File", sample(1:3, nrow(x), replace = TRUE))
colann <- data.frame(file = rep(paste0("File", 1:3), each = 10),
                     Channel = rep(names(x)[1:10], 3))

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
    ft2 <- readQFeatures(x, ecol = 1:10, name = "psms", fnames = 11)
    ft3 <- readQFeatures(x, colAnnotation = 1:10, name = "psms", fnames = 11)
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
    ft5 <- readSummarizedExperiment(x, ecol = 1:10, fnames = 11)
    expect_equal(ft3, ft5)
    ## Read data with only 1 quantitation column
    ft5 <- readSummarizedExperiment(x, ecol = 1, fnames = 11)
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


## ----------------------------------------------------------


test_that("readQFeatures: correct use", {
    #####################################
    ## Multiple batches
    qf <- readQFeatures(x, colann,
                        batchCol = "file",
                        channelCol = "Channel")
    expect_identical(sort(names(qf)), sort(unique(x$file)))
    expect_true(all(dims(qf)[2, ] == 10L))
    expect_true(sum(dims(qf)[1, ]) == nrow(x))
    ## Make sure all rownames start with "PSM"
    expect_true(all(grepl("^PSM", unlist(rownames(qf)))))
    ## Make sure the column names are as expected
    expectedCols <- paste0(rep(unique(x$file), 3),
                           rep(c("X126",  "X127C", "X127N", "X128C", "X128N",
                                 "X129C", "X129N", "X130C", "X130N", "X131" ), each = 3))
    expect_identical(character(0), setdiff(unlist(colnames(qf)), as.character(expectedCols)))
    #####################################
    ## Single batch
    onebatch <- x %>%
        dplyr::filter(file == "File1")
    qf <- readQFeatures(onebatch,
                        colann,
                        batchCol = "file",
                        channelCol = "Channel")
    expect_identical(dims(qf)[1, ],
                     c("File1" = nrow(onebatch)))
    expect_identical(dims(qf)[2, ],
                     c("File1" = 10L))
    #####################################
    ## Test remove empty columns
    qf <- readQFeatures(x, colann,
                        batchCol = "file",
                        channelCol = "Channel",
                        removeEmptyCols = TRUE)
    expect_identical(sort(names(qf)), sort(unique(x$file)))
    expect_true(all(dims(qf)[2, ] == rep(10, 3)))
    expect_true(sum(dims(qf)[1, ]) == nrow(x))
    #####################################
    ## Test suffix
    qf <- readQFeatures(x, colann,
                        batchCol = "file",
                        channelCol = "Channel",
                        suffix = paste0("_TMT", 1:10))
    expectedCols <- paste0(rep(unique(x$file), 10),
                           rep(paste0("_TMT", 1:10), each = 3))
    expectedCd <- DataFrame(colann)
    rownames(expectedCd) <- paste0(colann$file,
                                   "_TMT", 1:10)
    expect_true(all(unlist(colnames(qf)) %in% expectedCols))
    expect_identical(colData(qf)[sort(rownames(colData(qf))), ],
                     expectedCd[sort(rownames(expectedCd)), ])
    #####################################
    ## Test sep
    qf <- readQFeatures(x, colann,
                        batchCol = "file",
                        channelCol = "Channel",
                        suffix = paste0("TMT", 1:10),
                        sep = ".")
    expectedCols <- paste0(rep(unique(x$file), 10),
                           rep(paste0(".TMT", 1:10), each = 3))
    expectedCd <- DataFrame(colann)
    rownames(expectedCd) <- paste0(colann$file,
                                   ".TMT", 1:10)
    expect_true(all(unlist(colnames(qf)) %in% expectedCols))
    expect_identical(colData(qf)[sort(rownames(colData(qf))), ],
                     expectedCd[sort(rownames(expectedCd)), ])
})

test_that("readQFeatures: warnings", {
    ## Missing batch in metadata
    expect_warning(qf <- readQFeatures(x,
                                       dplyr::filter(colann,
                                                      file == "File1"),
                                        batchCol = "file",
                                        channelCol = "Channel"),
                   regexp = "Missing metadata. The features are removed")
    expect_identical(names(qf), "File1")
    expect_identical(dims(qf)[2, ], c("File1" = 10L))
    expect_identical(dims(qf)[1, ], c("File1" = sum(x$file == "File1")))
})

test_that("readQFeatures: error", {
    ## Suffix has not correct size
    expect_error(qf <- readQFeatures(x, colann,
                                     batchCol = "file",
                                     channelCol = "Channel",
                                     suffix = (1:2)),
                 regexp = "invalid rownames length")
    ## Suffix is not unique
    expect_error(expect_warning(
        qf <- readQFeatures(x, colann,
                            batchCol = "file",
                            channelCol = "Channel",
                            suffix = rep(1, 10)),
        regexp = "non-unique values"),
        regexp = "duplicate 'row.names'")
})
