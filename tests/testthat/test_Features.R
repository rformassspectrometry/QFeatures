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

test_that("updateObject", {
    data(feat3)
    ## Applying updateObject on feat3 should lead to the same object
    expect_identical(feat3, updateObject(feat3))
    ## Check verbose
    expect_message(updateObject(feat3, verbose = TRUE), 
                   regexp = "updateObject.*QFeatures")
})

test_that("addAssay: test correct use", {
    data(feat1)
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
    featS4 <- addAssay(feat1, el, dropColData = FALSE)
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

test_that("addAssay: test colData transfer", {
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
    colData(s3)$Group <- NULL
    ## Remove colData from assay (default)
    s3 <- addAssay(s3, s3[[1]], name = "assay3")
    expect_identical(colData(s3), colData(feat1))
    expect_true(isEmpty(colData(s3[["assay3"]])))
    ## Do not remove colData from assay (default)
    s3 <- addAssay(s3, s3[[1]], name = "assay3noColData", dropColData = FALSE)
    expect_identical(colData(s3), colData(feat1))
    expect_identical(colData(s3), colData(s3[["assay3noColData"]]))
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
    s7 <- addAssay(s7, se, name = "assay7")
    expect_identical(colData(s7), 
                     DataFrame(Group = 3:4,
                               row.names = c("S1", "S2")))
})
    

test_that("addAssay: test errors/warnings", {
    data(feat1)
    ## x is not a QFeatures
    expect_error(addAssay(feat1[[1]], y = feat1[[1]], name = "assay1"), 
                 regexp = "inherits.*QFeatures")
    ## y is corrupt
    corrupt <- feat1[[1]]
    corrupt@assays@data@listData[[1]] <- matrix()
    expect_error(addAssay(feat1, y = corrupt, name = "assay1"), 
                 regexp = "invalid.*SummarizedExperiment")
    ## name is ignored when y is provided as a list
    lse <- List(A = feat1[[1]], B = feat1[[1]])
    expect_warning(addAssay(feat1, y = lse, name = "foo"), 
                   regexp = "'name' is ignored")
    ## List of assays is unnamed
    ulse <- List(feat1[[1]], feat1[[1]])
    expect_error(addAssay(feat1, y = ulse), 
                 regexp = "named List")
    ## List must have unique names
    names(lse)[[2]] <- "A"
    expect_error(addAssay(feat1, y = lse), 
                 regexp = "names must be unique")
    ## Assay name already present in QFeatures
    names(lse)[[2]] <- "psms"
    expect_error(addAssay(feat1, y = lse), 
                 regexp = "already present")
    expect_error(addAssay(feat1, y = feat1[[1]], name = "psms"), 
                 regexp = "already present")
    ## One of the assays is not an SE
    lse <- List(A = feat1[[1]], B = matrix())
    expect_error(addAssay(feat1, y = lse), 
                 regexp = "inherit.*SummarizedExperiment")
    ## AssayLink is not associated to the assay
    expect_error(addAssay(feat1, y = feat1[[1]], name = "foo",
                          assayLinks = AssayLinks(names = "bar")), 
                 regexp = "assayLinks.*named after the assay")
})
# test_that("removeAssay", {
#     
# })
# 
# test_that("updateAssay", {
#     
# })
# 
# test_that("[[<-", {
#     
# })


test_that("[,QFeatures", {
    data(feat1)
    feat1 <- aggregateFeatures(feat1, "psms", "Sequence", "peptides")
    expect_true(expect_warning(validObject(feat1[, , "psms"]),
                               regexp = "'experiments' dropped; see 'metadata'"))
    expect_true(expect_warning(validObject(feat1[, , "peptides"]),
                               regexp = "'experiments' dropped; see 'metadata'"))
    expect_true(expect_warning(validObject(feat1[, , 1]),
                               regexp = "'experiments' dropped; see 'metadata'"))
    expect_true(expect_warning(validObject(feat1[, , 2]),
                               regexp = "'experiments' dropped; see 'metadata'"))
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
    hlpsms <- hlpsms[1:10, ]
    ft1 <- readQFeatures(hlpsms, ecol = 1:10, name = "psms", fname = "Sequence")
    rownames(ft1[[1]][1:2]) <- rep("1", 2)
    expect_error(validObject(ft1))
})


test_that("longFormat", {
    colData(feat2)$X <- 1:12
    ## Select a single colvars and rowvars
    lt <- longFormat(feat2, rowvars = "Prot", colvars = "X")
    ## Check dimensions
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+2L)
    ## Check content
    expect_identical(lt$Prot,
                     unname(unlist(lapply(rowData(feat2), 
                                          ## Repeat 4x because 4 samples
                                          function(x) rep(x$Prot, 4)))))
    ## Select a single colvars and no rowvars (make sure that 
    ## the implementation does not break the MAE implementation)
    lt <- longFormat(feat2, colvars = "X")
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+1L)
    ## Select multiple rowvars
    lt <- longFormat(feat2, rowvars = c("Prot", "x"))
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+2L)
    ## Test error: rowvars is missing in rowData
    expect_error(longFormat(feat2, rowvars = "y"),
                 regexp = "not found")
})


