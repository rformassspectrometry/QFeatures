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
                 regexp = "not compatible with dim.")
    ## Error: factor is too short
    expect_error(.splitSE(se, factor(1:3)),
                 regexp = "not compatible with dim")
})



test_that(".checkFilterVariables works", {
    rd <- List(X1 = DataFrame(x = 1:3,
                              y = 1:3),
               X2 = DataFrame(x = 1:3,
                              z = 1:3),
               X3 = DataFrame(x = 1:3,
                              z = 1:3))
    expect_true(all(.checkFilterVariables(rd, "x")))
    expect_message(.checkFilterVariables(rd, "x"),
                   "'x' found in 3 out of 3 assay")
    expect_equivalent(rowSums(.checkFilterVariables(rd, "z")), 2L)
    expect_message(.checkFilterVariables(rd, "z"),
                   "'z' found in 2 out of 3 assay")
    expect_equivalent(rowSums(.checkFilterVariables(rd, "y")), 1L)
    expect_message(.checkFilterVariables(rd, "y"),
                   "'y' found in 1 out of 3 assay")
    expect_error(.checkFilterVariables(rd, "v"),
                 "'v' is/are absent from all rowData.")
    w <- 1
    expect_error(.checkFilterVariables(rd, "w"),
                 "'w' is/are absent from all rowData.")
})
