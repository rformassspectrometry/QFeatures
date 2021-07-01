data("feat2")
feat3 <- joinAssays(feat2, i = 1:3, "joinedAssay")
feat3 <- aggregateFeatures(feat3, i = 4, fcol = "Prot", fun = colSums)


test_that("QFeatures-validity", {
    expect_true(validObject(feat3))
})

test_that(".valid_QFeatures_indices", {
    ## Rename assay to insert mismatch between assayLinks names and
    ## QFeatures names
    names(feat2@ExperimentList)[1] <- "foo"
    expect_error(validObject(feat2),
                 regexp = "Assay links names are wrong.")
})


test_that(".unique_row_names", {
    ## Rename assay rownanes to insert duplicated rownames
    rownames(feat2@ExperimentList[[1]])[1:3] <- rep("foo", 3)
    expect_error(validObject(feat2),
                 regexp = "duplicated row names")
})


test_that(".validAssayLink", {
    ## Valid assay link 1 parent
    expect_identical(.validAssayLink(feat3, 1), NULL)
    ## Valid assay link 3 parents
    expect_identical(.validAssayLink(feat3, 4), NULL)
    ## Valid assay link with empty parent
    test <- feat3[letters[1:3], , ]
    expect_identical(.validAssayLink(test, 3), NULL)
    
    ##-- Corrupt the @from slot --###
    ## Change @from slot in AssayLink to point to missing assay
    test@assayLinks[[1]]@from <- "foo"
    expect_error(.validAssayLink(test, 1),
                 regexp = "@from not valid")
    
    ##-- Corrupt the @hits slot with 1 parent (Hits object) --###
    ## Add hits that links a root assay from a missing assay
    test <- feat3
    test@assayLinks[[1]]@hits <- Hits(1:3, 1:3, 10, nrow(test[[1]]),
                                      names_from = letters[1:3],
                                      names_to = letters[1:3])
    expect_error(.validAssayLink(test, 1),
                 regexp = "point from a missing assay")
    ## Add hits that links an assay from missing features
    test <- feat3
    mcols(test@assayLinks[[5]]@hits)$names_from[[1]] <- "foo"
    expect_error(.validAssayLink(test, 5),
                 regexp = "point from missing features")
    ## Add hits that links an assay to missing features
    test <- feat3
    mcols(test@assayLinks[[5]]@hits)$names_to[[1]] <- "foo"
    expect_error(.validAssayLink(test, 5),
                 regexp = "point to missing features")
    
    ##-- Corrupt the @hits slot with multiple parents (List object) --###
    ## Add hits that links an assay to a missing features
    test <- feat3
    mcols(test@assayLinks[[4]]@hits[[2]])$names_from[[1]] <- "foo"
    expect_error(.validAssayLink(test, 4),
                 regexp = "point from missing features")
    ## Add hits that links an assay to a missing features
    test <- feat3
    mcols(test@assayLinks[[4]]@hits[[2]])$names_to[[1]] <- "foo"
    expect_error(.validAssayLink(test, 4),
                 regexp = "point to missing features")
})

test_that(".validAssayLinks", {
    ## Valid AssayLinks
    expect_identical(.validAssayLinks(feat3), NULL)
    ## Corrupt the AssayLinks @names by reverting the AssayLink order
    feat3@assayLinks <- feat3@assayLinks[rev(seq_along(feat3@assayLinks))]
    expect_error(.validAssayLinks(feat3),
                 "@names not valid")
})
