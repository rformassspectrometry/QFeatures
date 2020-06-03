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


test_that("subsetByFeatures_multiple", {
    ## Note this chunk tests the behavior of subsetting, but also the integrity
    ## of the  assaylinks (aggregation, one to one, multiple parents and 
    ## multiple child links) created during the processing
    se1 <- feat1[["psms"]][1:7, ]
    colnames(se1) <- paste0("Sample", 1:2)
    se2 <- feat1[["psms"]][3:10, ]
    colnames(se2) <- paste0("Sample", 3:4)
    fts <- Features(SimpleList(psms1 = se1,
                               psms2 = se2))
    fts <- joinAssays(fts, c("psms1", "psms2"), name = "psmsall")
    fts <- aggregateFeatures(fts, i = "psmsall", fcol = "Sequence",
                             name = "peptides")
    fts <- aggregateFeatures(fts, i = "peptides", fcol = "Protein",
                             name = "proteins")
    fts <- normalize(fts, i = "peptides", name = "normpeptides", method = "sum")
    fts <- aggregateFeatures(fts, i = "normpeptides", fcol = "Protein",
                             name = "normproteins")
    res1 <- subsetByFeature(fts, "ProtA")
    res2 <- fts["ProtA", ]
    expect_equal(res1, res2)
    expect_identical(dims(fts),
                     matrix(c(7L, 2L, 8L, 2L, 10L, 4L, 3L, 4L, 2L, 4L, 3L, 4L, 2L, 4L), 
                            nrow = 2, dimnames = list(NULL, names(fts))))
    expect_identical("ProtA",  
                     unique(unlist(lapply(experiments(res2), 
                                          function(x) rowData(x)$Protein))))
    
})
