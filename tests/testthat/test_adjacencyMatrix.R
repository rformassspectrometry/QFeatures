data(feat1)
feat1 <- aggregateFeatures(feat1, "psms", "Sequence", name = "peptides")
adj <- Matrix::sparseMatrix(c(1, 2, 3, 3),
                            c(1, 2, 1, 2), x = 1)
rownames(adj) <- rownames(feat1[[2]])
colnames(adj) <- paste0("Prot", c("A", "B"))

test_that("adjacencyMatrix,SummarizedExperiment setter and getter work", {
    se <- feat1[[2]]
    adjacencyMatrix(se) <- adj
    expect_identical(adjacencyMatrix(se), adj)
    expect_true("adjacencyMatrix" %in% names(rowData(se)))
    expect_error(adjacencyMatrix(se) <- adj)
    adjacencyMatrix(se, adjName = "adj2") <- adj
    expect_identical(adjacencyMatrix(se, adjName = "adj2"), adj)
    expect_true("adj2" %in% names(rowData(se)))
})

test_that("adjacencyMatrix,QFeatures setter and getter work", {
    adjacencyMatrix(feat1, 2) <- adj
    expect_identical(adjacencyMatrix(feat1[[2]]), adj)
    expect_identical(adjacencyMatrix(feat1, 2)[[1]], adj)
    expect_error(adjacencyMatrix(feat1, 2) <- adj)
    adjacencyMatrix(feat1, 2, adjName = "adj2") <- adj
    expect_identical(adjacencyMatrix(feat1[[2]], adjName = "adj2"), adj)
    expect_identical(adjacencyMatrix(feat1, 2, adjName = "adj2")[[1]], adj)
})
