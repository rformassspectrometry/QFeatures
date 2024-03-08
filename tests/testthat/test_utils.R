test_that(".splitSE", {
    m <- matrix(1:100, ncol = 10,
                dimnames = list(paste0("row", 1:10),
                                paste0("col", 1:10)))
    se <- SummarizedExperiment(assays = m,
                               rowData = DataFrame(rowDataCol = 1:nrow(m)%%3),
                               colData = DataFrame(colvar = 1:ncol(m)%%5))
    ## Split by row
    expect_identical(length(.splitSE(se, "rowDataCol")), 3L)
    ## Split by col
    expect_identical(length(.splitSE(se, "colvar")), 5L)
    ## Error: variable not found
    expect_error(.splitSE(se, "foo"),
                 regexp = "not found")
    ## Error: cannot split using more than 1 variable
    expect_error(.splitSE(se, c("Raw.file", "protein")),
                 regexp = "must be of lenght one")
    ## Error: factor is too short
    expect_error(.splitSE(se, factor(1:3)),
                 regexp = "not compatible with dim")
})
