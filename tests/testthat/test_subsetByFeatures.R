data(feat1)


test_that("subsetByFeatures", {
    feat1 <- aggregateFeatures(feat1, 1, fcol = "Sequence", name = "peptides", fun = colMedians)
    feat1 <- aggregateFeatures(feat1, 2, fcol = "Protein", name = "proteins", fun = colMedians)
    res1 <- subsetByFeature(feat1, "ProtA")
    res2 <- feat1["ProtA", ]
    expect_equal(res1, res2)
    expect_identical(lengths(feat1),
                     c(nrow(feat1[[1]]),
                       length(unique(rowData(feat1[[1]])[["Sequence"]])),
                       length(unique(rowData(feat1[[1]])[["Protein"]]))))
    expect_identical(sort(rownames(feat1[[2]])),
                     sort(unique(rowData(feat1[[1]])[["Sequence"]])))
    expect_identical(sort(rownames(feat1[[3]])),
                     sort(unique(rowData(feat1[[1]])[["Protein"]])))
})
