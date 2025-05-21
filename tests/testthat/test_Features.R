data(feat1)

test_that("empty QFeatures", {
    feat0 <- QFeatures()
    expect_true(validObject(feat0))
    expect_true(isEmpty(feat0))
    expect_null(show(feat0))
})

test_that("Manual QFeatures", {
    ## code from inst/scripts/test_data.R
    ## used to generate feat1
    psms <- matrix(1:20, ncol = 2)
    colnames(psms) <- paste0("S", 1:2)
    rowdata <- DataFrame(Sequence = c("SYGFNAAR", "SYGFNAAR", "SYGFNAAR", "ELGNDAYK",
                                      "ELGNDAYK", "ELGNDAYK", "IAEESNFPFIK",
                                      "IAEESNFPFIK", "IAEESNFPFIK", "IAEESNFPFIK"),
                         Protein = c("ProtA", "ProtA", "ProtA", "ProtA", "ProtA",
                                     "ProtA", "ProtB", "ProtB", "ProtB", "ProtB"),
                         Var = 1:10)
    rownames(rowdata) <- rownames(psms) <- paste0("PSM", 1:10)
    coldata <- DataFrame(Group = 1:2)
    rownames(coldata) <- colnames(psms)
    psms <- SummarizedExperiment(psms, rowData = rowdata)
    feat2 <- QFeatures(list(psms = psms), colData = coldata)
    expect_true(validObject(feat2))
    ## subsetting
    expect_null(show(feat2))
    expect_equal(psms, feat2[[1]])
    expect_equal(psms, feat2[["psms"]])
    expect_equal(feat2, feat2[1:10, 1:2, 1])
    ## compare to serialised data
    ## data(feat1)
    ## expect_equivalent(feat1, feat2)
})


test_that("Test indexing", {
    data(feat2)
    ## Test single indexing
    ## With numeric
    expect_identical(.normIndex(feat2, 1), "assay1")
    ## With character
    expect_identical(.normIndex(feat2, "assay1"), "assay1")
    ## With factor
    expect_identical(.normIndex(feat2, factor("assay1")), "assay1")
    ## With logical
    expect_identical(.normIndex(feat2, c(TRUE, FALSE, FALSE)), "assay1")

    ## Test multiple indexing
    ## With numeric
    expect_identical(.normIndex(feat2, 1:2), c("assay1", "assay2"))
    ## With character
    expect_identical(.normIndex(feat2, c("assay1", "assay2")),
                     c("assay1", "assay2"))
    ## With factor
    expect_identical(.normIndex(feat2, factor(c("assay1", "assay2"))),
                     c("assay1", "assay2"))
    ## With logical
    expect_identical(.normIndex(feat2, c(TRUE, TRUE, FALSE)),
                     c("assay1", "assay2"))

    ## Test wrong index
    ## Missing
    expect_error(.normIndex(feat2),
                 regexp = "missing, with no default")
    ## Logical is too short
    expect_error(.normIndex(feat2, TRUE),
                 regexp = "does not match the number")
    ## Index out of bounds
    expect_error(.normIndex(feat2, 4), regexp = "out of bounds")
    ## Assay name is absent...
    expect_error(.normIndex(feat2, "foo"),
                 regexp = "not found.*foo")
    ## ... unless requested
    expect_identical(.normIndex(feat2, "foo", allowAbsent = TRUE),
                     "foo")
})

test_that("updateObject", {
    data(feat3)
    ## Applying updateObject on feat3 should lead to the same object
    expect_identical(feat3, updateObject(feat3))
    ## Check verbose
    expect_message(updateObject(feat3, verbose = TRUE),
                   regexp = "updateObject.*QFeatures")
})

test_that("[[<-", {
    data("feat2")
    ## Check errors
    err <- feat2
    expect_error(err[[1:3]] <- experiments(feat3)[1:3],
                 regexp = "multiple replacement")
    expect_error(err[[1, j = 1]] <- err[[1]],
                 regexp = "invalid replacement")
    expect_error(err[[1, foo = 1]] <- err[[1]],
                 regexp = "invalid replacement")

    ## Check indexing
    charIndex <- numIndex <- feat2
    charIndex[["assay2"]] <- charIndex[["assay1"]]
    numIndex[[2]] <- numIndex[[1]]
    expect_identical(charIndex, numIndex)

    ## Scenario 1: use [[<- for replacement
    s1 <- feat2
    s1[[2]] <- s1[[1]]
    expect_identical(s1, replaceAssay(feat2, feat2[[1]], i = 2))

    ## Scenario 2: use [[<- for removal
    s2 <- feat2
    expect_warning(s2[[1]] <- NULL, regexp = "dropped")
    expect_identical(s2, expect_warning(removeAssay(feat2, i = 1),
                                        regexp = "dropped"))

    ## Scenario 3: use [[<- for adding
    s3 <- feat2
    s3[["foo"]] <- s3[[1]]
    expect_identical(s3, addAssay(feat2, feat2[[1]], name = "foo"))
})

test_that("dims,ncols,nrows", {
    ## Test dims
    data("feat3")
    expect_identical(dims(feat3),
                     matrix(c(7L, 8L, 10L, 3L, 2L, 3L, 2L,
                              rep(2L, 2), rep(4L, 5)),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(NULL, names(feat3))))
    ## Without use names
    expect_identical(dims(feat3, use.names = FALSE),
                     matrix(c(7L, 8L, 10L, 3L, 2L, 3L, 2L,
                              rep(2L, 2), rep(4L, 5)),
                            nrow = 2, byrow = TRUE))

    ## Test nrows
    expect_identical(dims(feat3)[1, ],
                     nrows(feat3))
    ## Without use names
    expect_identical(dims(feat3, use.names = FALSE)[1, ],
                     nrows(feat3, use.names = FALSE))

    ## Test ncols
    expect_identical(dims(feat3)[2, ],
                     ncols(feat3))
    ## Without use names
    expect_identical(dims(feat3, use.names = FALSE)[2, ],
                     ncols(feat3, use.names = FALSE))
})

test_that("coerce,MultiAssayExperiment,QFeatures-method", {
    data("feat3")
    mae <- as(feat3, "MultiAssayExperiment")
    expect_true(validObject(mae))
    qf <- as(mae, "QFeatures")
    exp <- feat3
    feat3@assayLinks <- AssayLinks(names = names(feat3))
    expect_identical(qf, feat3)
})

test_that("c,QFeatures-method", {
    data("feat3")

    ## Combine 2 QF objects
    suppressMessages(suppressWarnings(qf1 <- feat3[,, 1:2]))
    suppressMessages(suppressWarnings(qf2 <- feat3[,, 3:7]))
    expctd <- feat3
    expctd@assayLinks[[3]] <- AssayLink(names(feat3)[[3]]) ## this AssayLink will get lost
    expect_identical(c(qf1, qf2), expctd)

    ## Combine 3 QF objects
    suppressMessages(suppressWarnings(qf1 <- feat3[,, 1]))
    suppressMessages(suppressWarnings(qf2 <- feat3[,, 2]))
    suppressMessages(suppressWarnings(qf3 <- feat3[,, 3:7]))
    expect_identical(c(qf1, qf2, qf3), expctd)

    ## Error: Combine 1 QF object and 1 List object, or separate assays
    expect_error(c(qf1, experiments(feat3)),
                 regexp = "Consider using 'addAssay")
    expect_error(c(qf1, assay1 = feat3[[1]]),
                 regexp = "Consider using 'addAssay")
    ## Error: combine 1 QF object and 1 MAE object
    mae <- as(qf2, "MultiAssayExperiment")
    expect_error(c(qf1, mae), regexp = "with one or more MultiAssayExperiment")
    ## Error: combine 1 QF object and 1 object other than QF, MAE, SE or List
    expect_error(c(qf1, matrix()),
                 regexp = "coercing .matrix. to .QFeatures.")

    ## Warning: providing names is ignored
    expect_warning(cmbnd <- c(qf1, object1 = qf2, object2 = qf3),
                   regexp = "will be ignored")
    expect_identical(cmbnd, expctd)

    ## Check special cases with colData (up to now there was no colData)
    ## Same colData column in each object and same samples and same information
    qf3 <- qf1
    names(qf3) <- "psms3"
    qf1$foo <- qf3$foo <- "bar1"
    expect_identical(colData(c(qf1, qf3)),
                     DataFrame(foo = rep("bar1", 2),
                               row.names = paste0("Sample", 1:2)))
    ## Same colData column in each object and same samples and same information
    qf3$foo <- "bar2"
    expect_error(c(qf1, qf3), regexp = "conflict.*foo .in argument 2")
    ## Same colData column in each object and different samples
    qf2$foo <- "bar2"
    expect_identical(colData(c(qf1, qf2)),
                     DataFrame(foo = c("bar1", "bar1", "bar2", "bar2"),
                               row.names = paste0("Sample", 1:4)))
    ## Different colData column in each object and same samples
    qf3$foo <- NULL
    qf3$foo2 <- "bar2"
    expect_identical(colData(c(qf1, qf3)),
                     DataFrame(foo = rep("bar1", 2),
                               foo2 = rep("bar2", 2),
                               row.names = paste0("Sample", 1:2)))
    ## Different colData column in each object and different samples
    qf2$foo <- NULL
    qf2$foo2 <- "bar2"
    expect_identical(colData(c(qf1, qf2)),
                     DataFrame(foo = c("bar1", "bar1", NA, NA),
                               foo2 = c(NA, NA, "bar2", "bar2"),
                               row.names = paste0("Sample", 1:4)))
})

test_that("addAssay", {
    data(feat1)

    ## Check errors
    ## x is not a QFeatures
    expect_error(addAssay(feat1[[1]], y = feat1[[1]], name = "assay1"),
                 regexp = "inherits.*QFeatures")
    ## AssayLink is not associated to the assay
    expect_error(addAssay(feat1, y = feat1[[1]], name = "foo",
                          assayLinks = AssayLinks(names = "bar")),
                 regexp = "assayLinks.*named after the assay")
    ## Assay to add is corrupt
    corruptAssay <- feat1[[1]]
    corruptAssay@assays@data@listData[[1]] <- matrix(1)
    expect_error(addAssay(feat1, y = corruptAssay, name = "foo",
                          assayLinks = AssayLinks(names = "bar")),
                 regexp = "invalid class .SummarizedExperiment")

    ## Scenario 1: add an assay with same dimnames
    assay1 <- feat1[[1]]
    featS1 <- addAssay(feat1, y = assay1, name = "assay1")
    ## Check a new assay was added
    expect_identical(names(featS1), c("psms", "assay1"))
    ## Check the new assay contains the expected object
    expect_identical(featS1[[1]], featS1[[2]])
    ## Check the colData is unchanged
    expect_identical(colData(feat1), colData(featS1))
    ## Check the sampleMap is adapted correctly
    expect_identical(sampleMap(featS1),
                     rbind(sampleMap(feat1),
                           DataFrame(assay = rep("assay1", 2),
                                     primary = colnames(assay1),
                                     colname = colnames(assay1))))
    ## Check an AssayLinks object associated to the new assay is added
    expect_identical(names(featS1), names(featS1@assayLinks))

    ## Scenario 2: add an assay with a subset of the dimnames
    assay2 <- feat1[[1]][1:5, 1]
    featS2 <- addAssay(feat1, y = assay2, name = "assay2")
    ## Check a new assay was added
    expect_identical(names(featS2), c("psms", "assay2"))
    ## Check the new assay contains the expected object
    expect_identical(featS2[[1]][1:5, 1], featS2[[2]])
    ## Check the colData is unchanged
    expect_identical(colData(feat1), colData(featS2))
    ## Check the sampleMap is adapted correctly
    expect_identical(sampleMap(featS2),
                     rbind(sampleMap(feat1),
                           DataFrame(assay = "assay2",
                                     primary = colnames(assay2),
                                     colname = colnames(assay2))))
    ## Check an AssayLinks object associated to the new assay is added
    expect_identical(names(featS2), names(featS2@assayLinks))

    ## Scenario 3: add 2 assays with same dimnames
    assay3 <- feat1[[1]]
    el <- List(assayA = assay3, assayB = assay3)
    featS3 <- addAssay(feat1, el)
    ## Check the 2 assays were added and name is ignored
    expect_identical(names(featS3), c("psms", "assayA", "assayB"))
    ## Check the new assay contains the expected quantitative data.
    expect_identical(unname(assay(featS3[[1]])),
                     unname(assay(featS3[[2]])))
    ## Check the colData is adapted correctly
    expect_identical(colData(feat1), colData(featS3))
    ## Check the sampleMap is adapted correctly
    expect_identical(sampleMap(featS3),
                     rbind(sampleMap(feat1),
                           DataFrame(assay = rep(names(el), each = 2),
                                     primary = colnames(assay3),
                                     colname = colnames(assay3))))
    ## Check an AssayLinks object associated to the new assay is added
    expect_identical(names(featS2), names(featS2@assayLinks))

    ## Scenario 4: add 2 assays with different dimnames, one with colData
    ## the other without
    assay4 <- feat1[[1]]
    ## change row and sample names
    colnames(assay4) <- paste("foo", 1:ncol(assay4))
    rownames(assay4) <- paste("bar", 1:nrow(assay4))
    ## Add colData to one of the assays
    assay5 <- assay4
    colData(assay5)$foobar <- 1:2
    el <- List(assayA = assay4, assayB = assay5)
    featS4 <- addAssay(feat1, el)
    ## Check the colData is adapted correctly
    expect_identical(colData(featS4),
                     rbind(cbind(colData(feat1), foobar = NA),
                           DataFrame(foobar = 1:2, Group = NA,
                                     row.names = colnames(assay5))))
    ## Check the sampleMap is adapted correctly
    expect_identical(sampleMap(featS4),
                     rbind(sampleMap(feat1),
                           DataFrame(assay = rep(names(el), each = 2),
                                     primary = colnames(assay4),
                                     colname = colnames(assay4))))
    ## Test keeping the coldata
    featS4 <- addAssay(feat1, el)
    expect_true(!isEmpty(colData(featS4[["assayB"]])))

    ## Scenario 5: add an assay with assayLinks
    assay5 <- feat1[[1]]
    seq <- rowData(assay5)$Sequence
    al <- AssayLink("assay5", "psms", fcol = ".rows",
                    hits = findMatches(seq, seq))
    featS5 <- addAssay(feat1, y = assay5, name = "assay5",
                       assayLinks = al)
    ## Check the AssayLinks object associated to the new assay is added
    expect_identical(names(featS5), names(featS5@assayLinks))
    expect_identical(featS5@assayLinks[[2]], al)
})

test_that("replaceAssay", {
    data("feat2")

    ## Check errors
    ## x is not a QFeatures
    expect_error(replaceAssay(feat2[[1]], y = feat2[[1]], i = "psms1"),
                 regexp = "inherits.*QFeatures")
    ## Check indexing
    charIndex <- replaceAssay(feat2, feat2[[1]], i = "assay2")
    numIndex <- replaceAssay(feat2, feat2[[1]], i = 2)
    logIndex <- replaceAssay(feat2, feat2[[1]], i = names(feat2) == "assay2")
    expect_identical(charIndex, numIndex)
    expect_identical(charIndex, logIndex)

    ## Scenario 1: Replace an assay with itself should lead to an
    ## unmodified object
    expect_identical(feat1, replaceAssay(feat1, feat1[[1]], 1))
    ## But! when the colData in QFeatures is empty and the assays have
    ## non empty colData, then the replacement updates the colData
    s1 <- replaceAssay(feat2, experiments(feat2))
    expect_false(identical(s1, feat2))
    expect_identical(colData(s1),
                     rbind(colData(s1[[1]]),
                           colData(s1[[2]]),
                           cbind(colData(s1[[3]]), Var2 = NA)))
    ## Check the colData is still the same in each assay
    for (i in seq_along(experiments(s1))) {
        expect_identical(colData(s1[[i]]), colData(feat2[[i]]))
    }

    ## Scenario 2: Replace assay with colData and same samples.
    s2 <- replaceAssay(feat2, feat2[[1]], i = "assay2")
    expect_identical(s2[["assay2"]], feat2[[1]])

    ## Scenario 4: Replace assay that is parent and child of a single
    ## assay
    data("feat3")
    expect_warning(s4 <- replaceAssay(feat3, feat3[["psms1"]],
                                      i = "normpeptides"),
                   regexp = "Links between assays were lost")
    expect_identical(s4@assayLinks[["peptides"]],
                     feat3@assayLinks[["peptides"]])
    expect_identical(s4@assayLinks[["normpeptides"]],
                     AssayLink("normpeptides"))
    expect_identical(s4@assayLinks[["normproteins"]],
                     AssayLink("normproteins"))

    ## Scenario 5: Replace assay that is child of multiple assays
    expect_warning(s5 <- replaceAssay(feat3, feat3[["psms1"]],
                                      i = "psmsall"),
                   regexp = "Links between assays were lost")
    expect_identical(s5@assayLinks[["psmsall"]],
                     AssayLink("psmsall"))
    expect_identical(s5@assayLinks[["psms1"]],
                     AssayLink("psms1"))
    expect_identical(s5@assayLinks[["psms2"]],
                     AssayLink("psms2"))
    expect_identical(s5@assayLinks[["peptides"]],
                     AssayLink("peptides"))

    ## Scenario 6: Replace assay that is parent of multiple assays
    expect_warning(s6 <- replaceAssay(feat3, feat3[["psms1"]],
                                      i = "peptides"),
                   regexp = "Links between assays were lost")
    expect_identical(s6@assayLinks[["psmsall"]],
                     feat3@assayLinks[["psmsall"]])
    expect_identical(s6@assayLinks[["peptides"]],
                     AssayLink("peptides"))
    expect_identical(s6@assayLinks[["proteins"]],
                     AssayLink("proteins"))
    expect_identical(s6@assayLinks[["normpeptides"]],
                     AssayLink("normpeptides"))

    ## Scenario 7: Replace assay that is one of several parents, and
    ## is parent of no assays (hence no warning)
    expect_warning(s7 <- replaceAssay(feat3, feat3[["psms1"]],
                                      i = "psms2"),
                   regexp = "Links between assays were lost")
    expect_identical(s7@assayLinks[["psms2"]],
                     feat3@assayLinks[["psms2"]])
    expect_identical(s7@assayLinks[["psms1"]],
                     feat3@assayLinks[["psms1"]])
    expect_identical(s7@assayLinks[["psmsall"]],
                     QFeatures:::.create_assay_link(s7, from = "psms1",
                                                    to = "psmsall"))

    ## Scenario 8: multiple replacements, scenario 6 + 7
    el <- List(psms2 = feat3[["psms1"]],
               peptides = feat3[["psms1"]])
    expect_warning(s8 <- replaceAssay(feat3, el),
                   regexp = "Links between assays were lost")
    expect_identical(s8@assayLinks[["psms2"]],
                     AssayLink("psms2"))
    expect_identical(s8@assayLinks[["psmsall"]],
                     QFeatures:::.create_assay_link(s8, from = "psms1",
                                                    to = "psmsall"))
    expect_identical(s8@assayLinks[["peptides"]],
                     AssayLink("peptides"))
    expect_identical(s8@assayLinks[["proteins"]],
                     AssayLink("proteins"))
    expect_identical(s8@assayLinks[["normpeptides"]],
                     AssayLink("normpeptides"))

    ## Scenario 9: replace with a sample that has new column names
    se <- feat3[["psms1"]]
    colnames(se) <- cn <- paste0("foo", 1:ncol(se))
    expect_warning(s9 <- replaceAssay(feat3, se, i = "psms1"),
                   regexp = "Links between assays were lost")
    expect_identical(colnames(s9)[["psms1"]], cn)

    ## Scenario 10: replace with a sample that removes column names
    s10 <- replaceAssay(feat2, feat2[[2]], i = "assay1")
    expect_identical(rownames(colData(s10)), paste0("S", 5:12))
    expect_identical(unique(unlist(colnames(s10))), paste0("S", 5:12))

    ## Scenario 11: replacing a sample with the same dimnames doesn't
    ## remove feature links. More specifically, repacing an assay with
    ## itself should not change the QFeatures object
    expect_identical(replaceAssay(feat3, feat3[["psmsall"]], i = "psmsall"),
                     feat3)
})


test_that("removeAssay", {
    data("feat1")
    data("feat2")
    data("feat3")
    ## Check indexing
    expect_warning(charIndex <- removeAssay(feat2, i = "assay2"),
                   regexp = "dropped")
    expect_warning(numIndex <- removeAssay(feat2, i = 2),
                   regexp = "dropped")
    expect_warning(
        logIndex <- removeAssay(feat2, i = names(feat2) == "assay2"),
        regexp = "dropped")
    expect_identical(charIndex, numIndex)
    expect_identical(charIndex, logIndex)

    ## Scenario 1: remove the only assay
    expect_warning(s1 <- removeAssay(feat1, i = "psms"),
                   regexp = "dropped")
    expect_true(isEmpty(s1))

    ## Scenario 2: remove an assay that removes samples
    expect_warning(s2 <- removeAssay(feat2, i = "assay1"),
                   regexp = "dropped")
    expect_identical(unique(unlist(colnames(s2))), paste0("S", 5:12))
    expect_identical(rownames(colData(s2)), paste0("S", 5:12))

    ## Scenario 3: remove an assay that doesn't remove samples
    expect_warning(s3 <- removeAssay(feat3, i = "psms1"),
                   regexp = "dropped")
    expect_identical(colnames(s3),
                     colnames(feat3)[-1])
    expect_identical(rownames(colData(s3)),
                     paste0("Sample", 1:4))

    ## Scenario 4: remove multiple assays
    expect_warning(s2 <- removeAssay(feat2, i = 1:2),
                   regexp = "dropped")
    expect_identical(colnames(s2), CharacterList(assay3 = paste0("S", 9:12)))
    expect_identical(rownames(colData(s2)), paste0("S", 9:12))

    # Scenario 5: remove assay that is parent and child of a single assay
    expect_warning(s5 <- removeAssay(feat3, i = "normpeptides"),
                   regexp = "dropped")
    expect_false("normpeptides" %in% names(s5@assayLinks))
    expect_identical(s5@assayLinks[["normproteins"]],
                     AssayLink("normproteins"))

    ## Scenario 6: remove assay that is child of multiple assays
    expect_warning(s6 <- removeAssay(feat3, i = "psmsall"),
                   regexp = "dropped")
    expect_false("psmsall" %in% names(s6@assayLinks))
    expect_identical(s6@assayLinks[["psms1"]],
                     AssayLink("psms1"))
    expect_identical(s6@assayLinks[["psms2"]],
                     AssayLink("psms2"))
    expect_identical(s6@assayLinks[["peptides"]],
                     AssayLink("peptides"))

    ## Scenario 7: remove assay that is parent of multiple assays
    expect_warning(s7 <- removeAssay(feat3, i = "peptides"),
                   regexp = "dropped")
    expect_false("peptides" %in% names(s7@assayLinks))
    expect_identical(s7@assayLinks[["proteins"]],
                     AssayLink("proteins"))
    expect_identical(s7@assayLinks[["normpeptides"]],
                     AssayLink("normpeptides"))
})



test_that(".checkAssaysToInsert", {
    ## y is corrupt
    corrupt <- feat1[[1]]
    corrupt@assays@data@listData[[1]] <- matrix()
    expect_error(.checkAssaysToInsert(corrupt, feat1, name = "assay1"),
                 regexp = "invalid.*SummarizedExperiment")
    ## name is ignored when y is provided as a list
    lse <- List(A = feat1[[1]], B = feat1[[1]])
    expect_warning(.checkAssaysToInsert(lse, feat1, name = "foo"),
                   regexp = "'name' is ignored")
    ## List of assays is unnamed
    ulse <- List(feat1[[1]], feat1[[1]])
    expect_error(.checkAssaysToInsert(ulse, feat1),
                 regexp = "named List")
    ## List must have unique names
    names(lse)[[2]] <- "A"
    expect_error(.checkAssaysToInsert(lse, feat1),
                 regexp = "names must be unique")
    ## Assay name already present in QFeatures
    names(lse)[[2]] <- "psms"
    expect_error(.checkAssaysToInsert(lse, feat1),
                 regexp = "already present")
    expect_error(.checkAssaysToInsert(feat1[[1]], feat1, name = "psms"),
                 regexp = "already present")
    ## One of the assays is not an SE
    lse <- List(A = feat1[[1]], B = matrix())
    expect_error(.checkAssaysToInsert(lse, feat1),
                 regexp = "inherit.*SummarizedExperiment")

})

test_that("add/replaceAssay: test colData transfer", {
    data("feat1")
    ## Scenario 1: no colData in QFeatures, no colData in assay
    s1 <- feat1
    colData(s1)$Group <- NULL
    s1 <- addAssay(s1, s1[[1]], name = "assay1")
    expect_true(isEmpty(colData(s1)))
    ## Scenario 2: colData in QFeatures, no colData in assay
    s2 <- feat1
    s2 <- addAssay(s2, s2[[1]], name = "assay2")
    expect_identical(colData(s2), colData(feat1))
    ## Scenario 3: no colData in QFeatures, colData in assay
    s3 <- feat1
    colData(s3[[1]]) <- colData(s3)
    ## Do not remove colData from assay
    s3 <- addAssay(s3, s3[[1]], name = "assay3")
    expect_identical(colData(s3), colData(feat1))
    expect_identical(colData(s3), colData(s3[["assay3"]]))
    ## Scenario 4: colData in QFeatures, colData in assay with different
    ## samples and different colData variables
    s4 <- feat1
    se <- s4[[1]]
    colnames(se) <- paste0("foo", 1:ncol(se))
    se$bar <- letters[1:ncol(se)]
    s4 <- addAssay(s4, se, name = "assay4")
    expect_identical(colData(s4),
                     DataFrame(Group = c(1:2, NA, NA),
                               bar = c(NA, NA, "a", "b"),
                               row.names = c("S1", "S2", "foo1", "foo2")))
    ## Scenario 5: colData in QFeatures, colData in assay with same
    ## samples and different colData variables
    s5 <- feat1
    se <- s5[[1]]
    se$bar <- letters[1:ncol(se)]
    s5 <- addAssay(s5, se, name = "assay5")
    expect_identical(colData(s5),
                     DataFrame(Group = c(1:2),
                               bar = c("a", "b"),
                               row.names = c("S1", "S2")))
    ## Scenario 6: colData in QFeatures, colData in assay with different
    ## samples and same colData variables
    s6 <- feat1
    se <- s6[[1]]
    colnames(se) <- paste0("foo", 1:ncol(se))
    se$Group <- 1:2
    s6 <- addAssay(s6, se, name = "assay6")
    expect_identical(colData(s6),
                     DataFrame(Group = rep(1:2, 2),
                               row.names = c("S1", "S2", "foo1", "foo2")))
    ## Scenario 7: colData in QFeatures, colData in assay with same
    ## samples and same colData variables
    s7 <- feat1
    se <- s7[[1]]
    se$Group <- 3:4
    expect_error(addAssay(s7, se, name = "assay7"),
                 regexp = "colData in y have conflict.*Group.*assay7")
    ## Scenario 8: colData in QFeatures, no colData in replacement
    ## assay
    s8 <- feat1
    s8 <- replaceAssay(s8, s8[[1]], i = "psms")
    expect_identical(colData(s2), colData(feat1))
    ## Scenario 9: colData in QFeatures, colData in replacement assay
    ## same samples but different colData variables
    s9 <- feat1
    se <- s9[[1]]
    se$bar <- letters[1:ncol(se)]
    s9 <- replaceAssay(s9, se, i = "psms")
    expect_identical(colData(s9),
                     DataFrame(Group = c(1:2),
                               bar = c("a", "b"),
                               row.names = c("S1", "S2")))
    ## Scenario 10: colData in QFeatures, colData in replacement assay
    ## different samples and different colData variables
    s10 <- feat1
    se <- s10[[1]]
    se$bar <- letters[1:ncol(se)]
    colnames(se) <- paste0("foo", 1:ncol(se))
    s10 <- replaceAssay(s10, se, i = "psms")
    expect_identical(colData(s10),
                     DataFrame(Group = as.logical(c(NA, NA)),
                               bar = c("a", "b"),
                               row.names = c("foo1", "foo2")))
    ## Scenario 11: colData in QFeatures, colData in replacement assay
    ## different and common samples and different colData variables.
    ## Replacement adds new samples and removes old samples
    s11 <- feat1
    se <- s11[[1]]
    se$bar <- letters[1:ncol(se)]
    colnames(se)[[2]] <- "foo"
    s11 <- replaceAssay(s11, se, i = "psms")
    expect_identical(colData(s11),
                     DataFrame(Group = c(1L, NA),
                               bar = c("a", "b"),
                               row.names = c("S1", "foo")))
    ## Scenario 12: colData in QFeatures, colData in assay with same
    ## samples and same colData variables, but NA in QFeatures
    s12 <- feat1
    se <- s12[[1]]
    s12$Group <- NA
    se$Group <- 1:2
    s12 <- addAssay(s12, se, name = "assay7")
    expect_identical(colData(s12), colData(feat1))
})



test_that("[,QFeatures", {
    data(feat1)
    feat1 <- aggregateFeatures(feat1, "psms", "Sequence", "peptides")
    expect_true(expect_warning(validObject(feat1[, , "psms"]),
                               regexp = "'experiments' dropped; see 'drops"))
    expect_true(expect_warning(validObject(feat1[, , "peptides"]),
                               regexp = "'experiments' dropped; see 'drops"))
    expect_true(expect_warning(validObject(feat1[, , 1]),
                               regexp = "'experiments' dropped; see 'drops"))
    expect_true(expect_warning(validObject(feat1[, , 2]),
                               regexp = "'experiments' dropped; see 'drops"))
})


test_that("RowData", {
    data(feat2)
    rd <- rowData(feat2)
    expect_identical(names(rd), paste0("assay", 1:3))
    expect_identical(rd[[1]], rowData(feat2[[1]]))
    expect_identical(rd[[2]], rowData(feat2[[2]]))
    expect_identical(rd[[3]], rowData(feat2[[3]]))
})


test_that("RowData<-", {
    data(feat2)
    feat3 <- feat4 <- feat5 <- feat2
    value <- rowData(feat2)
    value[["assay1"]]$Prot <- letters[seq_len(nrow(value[["assay1"]]))]
    value[["assay1"]] <- value[["assay1"]][, "Prot",  drop = FALSE]
    value[["assay1"]]$foo <- rep("bar", nrow(value[["assay1"]]))
    rowData(feat3) <- value[-3]
    ## assay not in value are untouched
    expect_identical(feat3[[2]], feat2[[2]])
    ## replacing by untouched rowData leads to the same rowData
    expect_identical(feat3[[3]], feat2[[3]])
    ## rowvars in value and in rowData are replaced
    expect_identical(rowData(feat3[[1]])$Prot, letters[seq_len(nrow(feat3[[1]]))])
    ## rowvars in value but not in rowData are added
    expect_true("foo" %in% colnames(rowData(feat3[[1]])))
    expect_identical(rowData(feat3[[1]])$foo, rep("bar", nrow(feat3[[1]])))
    ## rowvars in rowData but not in value are untouched
    expect_identical(rowData(feat3[[1]])$x, rowData(feat2[[1]])$x)
    ## The value is not a DataFrameList
    value2 <- lapply(value, as.data.frame)
    rowData(feat4) <- value2
    expect_identical(feat3, feat4)
    ## invalide value leads to warning
    names(value) <- NULL
    expect_warning(rowData(feat5) <-  value,
                   regexp = "Could not find a common assay")
    expect_identical(feat5, feat2)
})


test_that("rowDataNames", {
    rdn <- rowDataNames(feat1)
    expect_identical(length(feat1), length(rdn))
    expect_identical(names(feat1), names(rdn))
    for (i in seq_along(length(feat1)))
        expect_identical(rdn[[i]], names(rowData(feat1[[i]])))
})


test_that("selectRowData", {
    x <- c("Sequence", "Protein")
    ft <- selectRowData(feat1, x)
    expect_identical(length(ft), length(feat1))
    expect_identical(names(ft), names(feat1))
    expect_identical(rowDataNames(ft)[[1]], x)
    expect_error(selectRowData(feat1))
    expect_message(ft <- selectRowData(feat1,
                                       c("Sequence", "Protein",
                                         "Var", "location", "pval",
                                         "var_not_found")))
    expect_identical(feat1, ft)
})


test_that("rbindRowData", {
    data(feat2)
    ## Rbind rowData from 1 assay
    rd <- rbindRowData(feat2, 1)
    expect_true(inherits(rd, "DFrame"))
    expect_identical(unique(rd$assay), names(feat2)[1])
    expect_identical(colnames(rd), c("assay", "rowname",
                                     colnames(rowData(feat2)[[1]])))
    expect_identical(nrow(rd), sum(dims(feat2)[1, 1]))
    ## Get all common variable from all assays
    rd <- rbindRowData(feat2, seq_along(feat2))
    expect_true(inherits(rd, "DFrame"))
    expect_identical(unique(rd$assay), names(feat2))
    expect_identical(colnames(rd), c("assay", "rowname", "Prot", "x"))
    expect_identical(nrow(rd), sum(dims(feat2)[1, ]))
    ## Warning no common variables
    a <- feat2[[1]]
    rowData(a) <- DataFrame(foo = "bar")
    feat3 <- addAssay(feat2, a, name = "assay4")
    expect_warning(rd <- rbindRowData(feat3, seq_along(feat3)),
                   regexp = "No common columns")
    expect_identical(DataFrame(), rd)
})


test_that("renaming", {
    data(feat1)
    feat1 <- aggregateFeatures(feat1, "psms", fcol = "Sequence",
                               name = "peptides", fun = colSums)
    feat1 <- aggregateFeatures(feat1, "peptides", fcol = "Protein",
                               name = "proteins", fun = colSums)
    expect_true(validObject(feat1))
    feat2 <- feat1
    names(feat2) <- LETTERS[1:3]
    expect_true(validObject(feat2))
    # Expect errors
    feat2@assayLinks[[1]]@name <- "foo"
    expect_error(validObject(feat2),
                 regexp = "@names not valid")
    feat2 <- feat1
    feat2@assayLinks[[1]]@from <- "bar"
    expect_error(validObject(feat2),
                 regexp = "@from not valid")
    expect_error(names(feat2) <- 1:3,
                 regexp = "must be a character")
    expect_error(names(feat2) <- letters[c(1,2,2)],
                 regexp = "is duplicated$")
    ## Test with assays that have multiple parents
    data("feat3")
    newnames <- paste0("foo", seq_along(feat3))
    names(feat3) <- newnames
    expect_true(validObject(feat3))
    expect_identical(names(feat3),  newnames)
})

## This does not seem to work when runing R CMD check...
# test_that("plotting", {
#     data(feat2)
#     feat2 <- joinAssays(feat2, i = 1:3)
#     feat2 <- aggregateFeatures(feat2, 4, "Prot", name = "proteins")
#     ## Plot QFeautres with interactive = FALSE
#     ## expect_doppelganger creates a snapshot and compares to a test
#     ## snapshot in /_snaps/
#     set.seed(1234)
#     vdiffr::expect_doppelganger("qFeatures-plot", plot(feat2, interactive = FALSE))
# })

test_that("assays must have unique rownames", {
    data(hlpsms)
    hlpsms <- hlpsms[1:10, ]
    ft1 <- readQFeatures(hlpsms, quantCol = 1:10, name = "psms", fname = "Sequence")
    ## Adapt in slots directly because our code doesn't allow anymore to run:
    # rownames(ft1[[1]][1:2]) <- rep("1", 2)
    rownames(ft1@ExperimentList@listData[[1]][1:2]) <- rep("1", 2)
    expect_error(validObject(ft1))
})


test_that("dropEmptyAssays", {
    ## Object is not of class QFeatures = error
    expect_error(
        dropEmptyAssays(matrix()),
        regexp = "QFeatures.*is not TRUE"
    )
    ## dims is out of bounds = error
    expect_error(
        dropEmptyAssays(QFeatures(), 1:3),
        regexp = "in '1:2'"
    )
    expect_error(
        dropEmptyAssays(QFeatures(), 3),
        regexp = "in '1:2'"
    )
    ## Test on empty QFeatures
    expect_identical(
        dropEmptyAssays(QFeatures(), 1:2),
        QFeatures()
    )
    ## Test on QFeatures with 1 non empty assays = no effecrt
    data("feat1")
    expect_identical(
        dropEmptyAssays(feat1, 1:2),
        feat1
    )
    ## Test on QFeatures with mulitple non empty assays (and assayLinks)
    data("feat3")
    expect_identical(
        dropEmptyAssays(feat3, 1:2),
        feat3
    )
    ## Test on QFeatures with 1 empty assay (no features)
    qf <- feat3
    expect_warning(
        qf[[1]] <- qf[[1]][numeric(), ],
        regexp = "Links between assays were lost"
    )
    test <- expect_warning(
        dropEmptyAssays(qf, 1),
        regexp = "'experiments' dropped"
    )
    exp <- suppressWarnings(qf[, , -1])
    expect_identical(test, exp)
    expect_identical(dropEmptyAssays(qf, 2), qf)
    ## Test on QFeatures with all empty assays (no features)
    qf <- feat1
    qf[[1]] <- qf[[1]][numeric(), ]
    test <- expect_warning(
        dropEmptyAssays(qf, 1),
        regexp = "'experiments' dropped"
    )
    expect_identical(test, QFeatures())
    expect_identical(dropEmptyAssays(qf, 2), qf)
    ## Test on QFeatures with 1 empty assay (no samples)
    data("feat3")
    qf <- feat3
    expect_warning(
        qf[[1]] <- qf[[1]][, numeric()],
        regexp = "Links between assays were lost"
    )
    expect_identical(dropEmptyAssays(qf, 1), qf)
    test <- expect_warning(
        dropEmptyAssays(qf, 2),
        regexp = "'experiments' dropped"
    )
    exp <- suppressWarnings(qf[, , -1])
    expect_identical(test, exp)
    ## Test on QFeatures with all empty assays (no samples)
    qf <- feat1
    qf[[1]] <- qf[[1]][, numeric()]
    expect_identical(dropEmptyAssays(qf, 1), qf)
    test <- expect_warning(
        dropEmptyAssays(qf, 2),
        regexp = "'experiments' dropped"
    )
    expect_identical(test, QFeatures())
})

test_that("setQFeaturesType", {
    qf <- QFeatures()
    expect_equal(metadata(qf)[["._type"]], NULL)
    qf <- setQFeaturesType(qf, "bulk")
    expect_equal(metadata(qf)[["._type"]], "bulk")
    qf <- setQFeaturesType(qf, "scp")
    expect_equal(metadata(qf)[["._type"]], "scp")
    expect_error(setQFeaturesType(qf, "invalid"), "Invalid QFeatures type")
})

test_that("getQFeaturesType", {
    qf <- QFeatures()
    expect_message(type <- getQFeaturesType(qf),
                   "No explicit type set for this QFeatures object")
    expect_equal(type, "bulk")
    qf <- setQFeaturesType(qf, "scp")
    expect_equal(getQFeaturesType(qf), "scp")
})

test_that("validQFeaturesTypes", {
    expect_equal(validQFeaturesTypes(), c("bulk", "scp"))
})
