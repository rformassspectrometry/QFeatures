data(feat2)

.sampleAggregationSE <- function(group = c("a", "a", "b", "b")) {
    se <- feat2[["assay1"]]
    colData(se)$group <- group
    colData(se)$batch <- c("x", "y", "x", "y")
    colData(se)$fixed <- rep(c("A", "B"), each = 2)
    colData(se)$variable <- c("u", "v", "u", "v")
    rowData(se)$feature <- rownames(se)
    se
}

.expectedByGroup <- function(se) {
    m <- assay(se)
    cbind(a = rowMeans(m[, 1:2]),
          b = rowMeans(m[, 3:4]))
}

.expectedByGroupBatch <- function(se) {
    ans <- assay(se)
    storage.mode(ans) <- "double"
    colnames(ans) <- paste(colData(se)$group, colData(se)$batch, sep = "_")
    ans
}

test_that("aggregateSamples dispatches for SummarizedExperiment and QFeatures", {
    se <- .sampleAggregationSE()
    expected <- .expectedByGroup(se)

    se2 <- aggregateSamples(se, "group", fun = rowMeans, moreFun = list())

    expect_s4_class(se2, "SummarizedExperiment")
    expect_identical(assay(se2), expected)

    qf <- QFeatures(list(raw = se))
    qf2 <- aggregateSamples(qf, "raw", "group",
                            name = "grouped", fun = rowMeans,
                            moreFun = list())

    expect_s4_class(qf2, "QFeatures")
    expect_named(qf2, c("raw", "grouped"))
    expect_identical(assay(qf2[["grouped"]]), expected)
})

test_that("aggregateSamples,SummarizedExperiment: metadata and extra assays", {
    se <- .sampleAggregationSE()

    se2 <- aggregateSamples(se, "group", fun = rowMeans)

    expect_named(assays(se2), c("assay", "aggcounts"))
    expect_equal(assay(se2, "aggcounts"),
                 matrix(2, nrow = nrow(se), ncol = 2,
                        dimnames = list(rownames(se), c("a", "b"))))
    expect_identical(rowData(se2), rowData(se))
    expect_named(colData(se2), c("group", "fixed"))
    expect_identical(colData(se2)$group, c("a", "b"))
    expect_identical(colData(se2)$fixed, c("A", "B"))
    expect_true(validObject(se2))

    se3 <- aggregateSamples(se, "group", fun = rowMeans,
                            moreFun = list(aggsums = rowSums))

    expect_named(assays(se3), c("assay", "aggsums"))
    expect_equal(assay(se3, "aggsums"),
                 cbind(a = rowSums(assay(se)[, 1:2]),
                       b = rowSums(assay(se)[, 3:4])))
})

test_that("aggregateSamples,SummarizedExperiment: multiple scol columns", {
    se <- .sampleAggregationSE()
    se2 <- aggregateSamples(se, c("group", "batch"),
                            fun = rowMeans, moreFun = list())

    expect_s4_class(se2, "SummarizedExperiment")
    expect_identical(assay(se2), .expectedByGroupBatch(se))
})

test_that("aggregateSamples,SummarizedExperiment: unused factor levels", {
    se <- .sampleAggregationSE(
        group = factor(c("a", "a", "b", "b"), levels = c("a", "b", "c"))
    )
    se2 <- aggregateSamples(se, "group", fun = rowMeans, moreFun = list())

    expect_identical(colnames(se2), c("a", "b"))
    expect_identical(rownames(colData(se2)), c("a", "b"))
    expect_equal(assay(se2), .expectedByGroup(se))
    expect_true(validObject(se2))
})

test_that("aggregateSamples,SummarizedExperiment: errors", {
    se <- .sampleAggregationSE()

    expect_error(aggregateSamples(se, fun = rowMeans, moreFun = list()),
                 "'scol' is required.", fixed = TRUE)
    expect_error(aggregateSamples(se, character(), fun = rowMeans,
                                  moreFun = list()),
                 "'scol' is required.", fixed = TRUE)
    expect_error(aggregateSamples(se, 1, fun = rowMeans, moreFun = list()),
                 "'scol' must be a non-missing character vector.",
                 fixed = TRUE)
    expect_error(aggregateSamples(se, NA_character_, fun = rowMeans,
                                  moreFun = list()),
                 "'scol' must be a non-missing character vector.",
                 fixed = TRUE)
    expect_error(aggregateSamples(se, "missing", fun = rowMeans,
                                  moreFun = list()),
                 "'scol' not found in the assay's colData.",
                 fixed = TRUE)
    expect_error(aggregateSamples(se, c("group", "missing"),
                                  fun = rowMeans, moreFun = list()),
                 "'scol' not found in the assay's colData.",
                 fixed = TRUE)
})

test_that("aggregateSamples,QFeatures: argument errors", {
    se <- .sampleAggregationSE()
    qf <- QFeatures(list(raw = se))

    expect_error(aggregateSamples(qf, scol = "group"),
                 "'i' must be a non-missing character vector of length 1.",
                 fixed = TRUE)
    expect_error(aggregateSamples(qf, 1, "group"),
                 "'i' must be a non-missing character vector of length 1.",
                 fixed = TRUE)
    expect_error(aggregateSamples(qf, c("raw", "raw"), "group"),
                 "'i' must be a non-missing character vector of length 1.",
                 fixed = TRUE)
    expect_error(aggregateSamples(qf, "raw", "group", name = 1),
                 "'name' must be a non-missing character vector of length 1.",
                 fixed = TRUE)
    expect_error(aggregateSamples(qf, "raw", "group", name = c("one", "two")),
                 "'name' must be a non-missing character vector of length 1.",
                 fixed = TRUE)
    expect_error(aggregateSamples(qf, "raw", "group", name = "raw"),
                 "There's already an assay named: 'raw'.",
                 fixed = TRUE)
    expect_error(aggregateSamples(qf, "missing", "group", name = "grouped"),
                 "The following assay(s) is/are not found:missing",
                 fixed = TRUE)
})

test_that("aggregateSamples,QFeatures: multiple scol columns", {
    se <- .sampleAggregationSE()
    qf <- QFeatures(list(raw = se))
    qf2 <- aggregateSamples(qf, "raw", c("group", "batch"),
                            name = "grouped", fun = rowMeans,
                            moreFun = list())

    expect_s4_class(qf2, "QFeatures")
    expect_identical(assay(qf2[["grouped"]]), .expectedByGroupBatch(se))
    expect_true(validObject(qf2))
})

test_that("aggregateSamples,QFeatures: arguments to fun", {
    se <- .sampleAggregationSE()
    assay(se)[1, 1] <- NA
    qf <- QFeatures(list(raw = se))

    se2 <- aggregateSamples(se, "group", fun = rowMeans,
                            moreFun = list(), na.rm = TRUE)
    qf2 <- aggregateSamples(qf, "raw", "group", name = "grouped",
                            fun = rowMeans, moreFun = list(),
                            na.rm = TRUE)

    expect_equal(assay(qf2[["grouped"]]), assay(se2))
    expect_true(validObject(qf2))
})

test_that("aggregateSamples,QFeatures: scol errors", {
    se <- .sampleAggregationSE(group = c("S1", "S1", "S2", "S2"))
    qf <- QFeatures(list(raw = se))

    expect_error(aggregateSamples(qf, "raw"),
                 "'scol' is required.", fixed = TRUE)
    expect_error(aggregateSamples(qf, "raw", character(), name = "grouped",
                                  fun = rowMeans, moreFun = list()),
                 "'scol' is required.", fixed = TRUE)
    expect_error(aggregateSamples(qf, "raw", 1, name = "grouped",
                                  fun = rowMeans, moreFun = list()),
                 "'scol' must be a non-missing character vector.",
                 fixed = TRUE)
    expect_error(aggregateSamples(qf, "raw", NA_character_, name = "grouped",
                                  fun = rowMeans, moreFun = list()),
                 "'scol' must be a non-missing character vector.",
                 fixed = TRUE)
})

test_that("aggregateSamples,QFeatures: addAssay handles colData rows appending", {
    se <- .sampleAggregationSE(group = c("S1", "S1", "S2", "S2"))
    qf <- QFeatures(list(raw = se))

    qf2 <- aggregateSamples(qf, "raw", "group", name = "grouped",
                            fun = rowMeans, moreFun = list())

    expect_named(qf2, c("raw", "grouped"))
    expect_identical(colnames(qf2[["grouped"]]), c("S1", "S2"))
    expect_true(all(c("S1", "S2") %in% rownames(colData(qf2))))
    expect_true(validObject(qf2))
})

test_that("aggregateSamples,QFeatures: reuses compatible colData rows", {
    se <- .sampleAggregationSE()
    qf <- QFeatures(list(raw1 = se, raw2 = se))

    qf2 <- suppressWarnings(
        aggregateSamples(qf, "raw1", "group", name = "grouped1",
                         fun = rowMeans, moreFun = list())
    )
    qf2 <- suppressWarnings(
        aggregateSamples(qf2, "raw2", "group", name = "grouped2",
                         fun = rowMeans, moreFun = list())
    )

    expect_named(qf2, c("raw1", "raw2", "grouped1", "grouped2"))
    expect_equal(assay(qf2[["grouped1"]]), .expectedByGroup(se))
    expect_equal(assay(qf2[["grouped2"]]), .expectedByGroup(se))
    expect_true(validObject(qf2))
})

test_that("aggregateSamples,QFeatures: assays and links", {
    se <- .sampleAggregationSE()
    qf <- QFeatures(list(raw = se))
    qf2 <- aggregateSamples(qf, "raw", "group", name = "grouped",
                            fun = rowMeans,
                            moreFun = list(aggsums = rowSums))
    alink <- assayLink(qf2, "grouped")

    expect_named(assays(qf2[["grouped"]]), c("assay", "aggsums"))
    expect_equal(assay(qf2[["grouped"]], "aggsums"),
                 cbind(a = rowSums(assay(se)[, 1:2]),
                       b = rowSums(assay(se)[, 3:4])))
    expect_identical(rowData(qf2[["grouped"]]), rowData(se))
    expect_identical(alink@name, "grouped")
    expect_identical(alink@from, "raw")
    expect_identical(alink@fcol, "._oneToOne")
    expect_identical(alink@hits@from, seq_len(nrow(se)))
    expect_identical(alink@hits@to, seq_len(nrow(se)))
    expect_true(validObject(qf2))

    qf3 <- aggregateSamples(qf, "raw", "group", name = "counted",
                            fun = rowMeans)

    expect_named(assays(qf3[["counted"]]), c("assay", "aggcounts"))
    expect_equal(assay(qf3[["counted"]], "aggcounts"),
                 matrix(2, nrow = nrow(se), ncol = 2,
                        dimnames = list(rownames(se), c("a", "b"))))
    expect_true(validObject(qf3))
})
