test_that("aggregateSamples dispatches for SummarizedExperiment and QFeatures", {
    m <- matrix(seq_len(8),
                nrow = 2,
                dimnames = list(c("f1", "f2"), paste0("s", seq_len(4))))
    se <- SummarizedExperiment(
        assays = list(assay = m),
        colData = DataFrame(group = c("a", "a", "b", "b"))
    )
    expected <- cbind(a = rowMeans(m[, 1:2]),
                      b = rowMeans(m[, 3:4]))

    se2 <- aggregateSamples(se, "group", fun = rowMeans, moreFun = list())

    expect_s4_class(se2, "SummarizedExperiment")
    expect_identical(assay(se2), expected)

    qf <- QFeatures(list(raw = se))
    invisible(capture.output(
        qf2 <- aggregateSamples(qf, "raw", "group",
                                name = "grouped", fun = rowMeans,
                                moreFun = list())
    ))

    expect_s4_class(qf2, "QFeatures")
    expect_named(qf2, c("raw", "grouped"))
    expect_identical(assay(qf2[["grouped"]]), expected)
})

test_that("aggregateSamples validates QFeatures i and name", {
    m <- matrix(seq_len(8),
                nrow = 2,
                dimnames = list(c("f1", "f2"), paste0("s", seq_len(4))))
    se <- SummarizedExperiment(
        assays = list(assay = m),
        colData = DataFrame(group = c("a", "a", "b", "b"))
    )
    qf <- QFeatures(list(raw = se))

    expect_error(
        aggregateSamples(qf, scol = "group"),
        "'i' must be a non-missing character vector of length 1.",
        fixed = TRUE
    )
    expect_error(
        aggregateSamples(qf, 1, "group"),
        "'i' must be a non-missing character vector of length 1.",
        fixed = TRUE
    )
    expect_error(
        aggregateSamples(qf, c("raw", "raw"), "group"),
        "'i' must be a non-missing character vector of length 1.",
        fixed = TRUE
    )
    expect_error(
        aggregateSamples(qf, "raw", "group", name = 1),
        "'name' must be a non-missing character vector of length 1.",
        fixed = TRUE
    )
    expect_error(
        aggregateSamples(qf, "raw", "group", name = c("one", "two")),
        "'name' must be a non-missing character vector of length 1.",
        fixed = TRUE
    )
})
