data(feat3)

test_that("subsetByFeatures", {
    feat1 <- aggregateFeatures(feat1, 1, fcol = "Sequence", name = "peptides", fun = colMedians)
    feat1 <- aggregateFeatures(feat1, 2, fcol = "Protein", name = "proteins", fun = colMedians)
    ## Subset a feature that is present in all assays
    res1 <- subsetByFeature(feat1, "ProtA")
    res2 <- feat1["ProtA", ]
    res3 <- feat1["ProtA", , ]
    expect_identical(res1, res2)
    expect_identical(res1, res3)
    expect_identical(lengths(feat1),
                     c(nrow(feat1[[1]]),
                       length(unique(rowData(feat1[[1]])[["Sequence"]])),
                       length(unique(rowData(feat1[[1]])[["Protein"]]))))
    expect_identical(sort(rownames(feat1[[2]])),
                     sort(unique(rowData(feat1[[1]])[["Sequence"]])))
    expect_identical(sort(rownames(feat1[[3]])),
                     sort(unique(rowData(feat1[[1]])[["Protein"]])))
    ## Subset a feature that is present in some assays
    res1 <- subsetByFeature(feat1, "SYGFNAAR")
    res2 <- feat1["SYGFNAAR", ]
    res3 <- feat1["SYGFNAAR", , ]
    expect_identical(res1, res2)
    expect_identical(res1, res3)
    ## Check the QFeatures subsetting follows the same behavior as
    ## MultiAssayExperiment
    expect_identical(dims(subsetByFeature(feat1, "PSM1")),
                     dims(subsetByRow(feat1, "PSM1")))
})

test_that("subsetByFeatures: full pipeline", {
    ## Subsetting "ProtA" will go through all assays as they all contains the
    ## "Protein" variable in the `rowData`
    ftsub <- subsetByFeature(feat3, "ProtA")
    expect_true(validObject(ftsub))
    expect_identical(ftsub, feat3["ProtA", ])
    expect_identical(dims(ftsub),
                     matrix(c(6L, 2L, 4L, 2L, 6L, 4L, 2L, 4L, 1L, 4L, 2L, 4L, 1L, 4L),
                            nrow = 2,
                            dimnames = list(NULL, c("psms1", "psms2", "psmsall", "peptides", "proteins", "normpeptides", "normproteins"))))
    expect_identical("ProtA",
                     unique(unlist(lapply(experiments(ftsub),
                                          function(x) rowData(x)$Protein))))

    ## Subsetting "SYGFNAAR" will subset only the peptides and psms assays as
    ## the protein assays do not contain the peptide "Sequence" variable
    ftsub <- feat3["SYGFNAAR", ]
    expect_identical(length(ftsub), 7L)
    expect_identical(names(ftsub), names(feat3))
    expect_identical(dims(ftsub)[2, ], dims(feat3)[2, ])
    expect_identical(dims(ftsub)[1, ], dims(feat3)[1, ] - c(4L, 7L, 7L, rep(2L, 4)))
    expect_identical("SYGFNAAR",
                     unique(unlist(lapply(rowData(ftsub),
                                          function(x) x$Sequence))))
})
