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
    expect_true(all(.checkFilterVariables(rd, "x", keep = FALSE)))
    expect_true(all(.checkFilterVariables(rd, "x", keep = TRUE)))
    expect_message(.checkFilterVariables(rd, "x", keep = FALSE),
                   "'x' found in 3 out of 3 assay")
    expect_message(.checkFilterVariables(rd, "x", keep = TRUE),
                   "'x' found in 3 out of 3 assay")
    expect_equivalent(rowSums(.checkFilterVariables(rd, "z", keep = TRUE)), 2L)
    expect_equivalent(rowSums(.checkFilterVariables(rd, "z", keep = FALSE)), 2L)
    expect_message(.checkFilterVariables(rd, "z", keep = FALSE),
                   "'z' found in 2 out of 3 assay")
        expect_message(.checkFilterVariables(rd, "z", keep = TRUE),
                   "'z' found in 2 out of 3 assay")
    expect_equivalent(rowSums(.checkFilterVariables(rd, "y", keep = FALSE)), 1L)
    expect_equivalent(rowSums(.checkFilterVariables(rd, "y", keep = TRUE)), 1L)
    expect_message(.checkFilterVariables(rd, "y", keep = FALSE),
                   "'y' found in 1 out of 3 assay")
    expect_message(.checkFilterVariables(rd, "y", keep = TRUE),
                   "'y' found in 1 out of 3 assay")
    expect_error(.checkFilterVariables(rd, "v", keep = FALSE),
                 "'v' is/are absent from all rowData.")
    expect_error(.checkFilterVariables(rd, "v", keep = TRUE),
                 "'v' is/are absent from all rowData.")
    w <- 1
    expect_error(.checkFilterVariables(rd, "ww", keep = FALSE),
                 "'ww' is/are absent from all rowData.")
    expect_error(.checkFilterVariables(rd, "ww", keep = TRUE),
                 "'ww' is/are absent from all rowData.")
})


test_that(".setAssayRownames works", {
    data(feat2)
    newnames <- CharacterList(LETTERS[1:10],
                              LETTERS[1:4],
                              LETTERS[1:7])
    for (i in 1:3)
        rowData(feat2[[i]])[["NEWNAME"]] <- newnames[[i]]
    ans <- .setAssayRownames(feat2, "NEWNAME")
    expect_true(validObject(ans))
    expect_identical(unname(rownames(ans)),
                     newnames)
})

test_that(".setAssayRownames works with duplicated names", {
    data(feat2)
    newnames <- CharacterList(c(LETTERS[1:5], LETTERS[1:5]),
                              LETTERS[1:4],
                              LETTERS[1:7])
    for (i in 1:3)
        rowData(feat2[[i]])[["NEWNAME"]] <- newnames[[i]]
    ans <- .setAssayRownames(feat2, "NEWNAME")
    expect_true(validObject(ans))
    ## all but first names are the same
    expect_identical(unname(rownames(ans)[-1]),
                     newnames[-1])
    expect_identical(unname(rownames(ans)[[1]]),
                     make.unique(newnames[[1]]))
})