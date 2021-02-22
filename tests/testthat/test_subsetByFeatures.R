## Note the dataset is a nice example that tests QFeatures as it could routinely
## be used within a pipeline. Note that the 3 special cases of AssayLinks are
## present in this dataset:
##     * One to one link: link between 1 parent and 1 child with one to one row
##       mapping. E.g. "peptides" to "normpeptides"
##     * One parent to multiple children: link between 1 parents and multiple
##       children mapping. E.g. "peptides" to "normpeptides" and "proteins"
##     * Multiple parents to one child: link between multiple parents and one
##       child mapping. E.g. "psms1" and "psms2" to "psmsall"
##
## psms1 ---
##           \
##             ----> psmsall ----> peptides ----> proteins
##           /                        |
## psms2 ---                          |
##                               normpeptides ----> normproteins
##
## of the  assaylinks (aggregation, one to one, multiple parents and
## multiple child links) created during the processing
data(feat1)
se1 <- feat1[["psms"]][1:7, ]
colnames(se1) <- paste0("Sample", 1:2)
se2 <- feat1[["psms"]][3:10, ]
colnames(se2) <- paste0("Sample", 3:4)
fts <- QFeatures(SimpleList(psms1 = se1,
                           psms2 = se2))
fts <- joinAssays(fts, c("psms1", "psms2"), name = "psmsall")
fts <- aggregateFeatures(fts, i = "psmsall", fcol = "Sequence",
                         name = "peptides")
fts <- aggregateFeatures(fts, i = "peptides", fcol = "Protein",
                         name = "proteins")
fts <- normalize(fts, i = "peptides", name = "normpeptides", method = "sum")
fts <- aggregateFeatures(fts, i = "normpeptides", fcol = "Protein",
                         name = "normproteins")

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
    res1 <- expect_warning(subsetByFeature(feat1, "SYGFNAAR"),
                           regexp = "'experiments' dropped; see 'metadata'")
    res2 <- expect_warning(feat1["SYGFNAAR", ],
                           regexp = "'experiments' dropped; see 'metadata'")
    res3 <- expect_warning(feat1["SYGFNAAR", , ],
                           regexp = "'experiments' dropped; see 'metadata'")
    expect_identical(res1, res2)
    expect_identical(res1, res3)
})

test_that("subsetByFeatures: full pipeline", {
    ## Subsetting "ProtA" will go through all assays as they all contains the
    ## "Protein" variable in the `rowData`
    ftsub <- subsetByFeature(fts, "ProtA")
    expect_identical(ftsub, fts["ProtA", ])
    expect_identical(dims(ftsub),
                     matrix(c(6L, 2L, 4L, 2L, 6L, 4L, 2L, 4L, 1L, 4L, 2L, 4L, 1L, 4L),
                            nrow = 2,
                            dimnames = list(NULL, c("psms1", "psms2", "psmsall", "peptides", "proteins", "normpeptides", "normproteins"))))
    expect_identical("ProtA",
                     unique(unlist(lapply(experiments(ftsub),
                                          function(x) rowData(x)$Protein))))

    ## Subsetting "SYGFNAAR" will subset only the peptides and psms assays as
    ## the protein assays do not contain the peptide "Sequence" variable
    expect_warning(expect_message(ftsub <- fts["SYGFNAAR", ],
                                  regexp = "removing 8 sampleMap rows not in names"),
                   regexp = "'experiments' dropped; see 'metadata'")
    expect_identical(length(ftsub), 5L)
    expect_identical(dims(ftsub),
                     matrix(c(3L, 2L, 1L, 2L, 3L, 4L, 1L, 4L, 1L, 4L),
                            nrow = 2,
                            dimnames = list(NULL, c("psms1", "psms2", "psmsall", "peptides", "normpeptides"))))
    expect_identical("SYGFNAAR",
                     unique(unlist(lapply(experiments(ftsub),
                                          function(x) rowData(x)$Sequence))))
})
