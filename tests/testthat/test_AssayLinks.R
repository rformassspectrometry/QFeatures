data(feat1)
feat1 <- aggregateFeatures(feat1, "psms", "Sequence", "peptides")
featMultiP <- addAssay(feat1, feat1[["psms"]], name = "psms2")
alMulti <- .create_assay_link_multi(featMultiP, 
                                    from = c("psms", "psms2"),
                                    to = "peptides",  
                                    varsFrom = c("Sequence", "Sequence"),
                                    varTo = "Sequence")

test_that(".create_assay_link", {
    ## Test correct usage
    ## we here use expect_equivalent and not identical  because the `hits` in 
    ## `feat1` are sorted by `Sequence` by `aggregateFeatures`
    expect_equivalent(feat1@assayLinks[[2]], 
                      .create_assay_link(rdFrom = rowData(feat1[["psms"]]), 
                                         rdTo = rowData(feat1[["peptides"]]),
                                         from = "psms", to = "peptide",
                                         varFrom = "Sequence",                                          
                                         varTo = "Sequence"))
    ## Test errors
    ## No match between feature variables
    expect_error(.create_assay_link(rdFrom = rowData(feat1[["psms"]]), 
                                    rdTo = rowData(feat1[["peptides"]]),
                                    from = "psms", to = "peptide",
                                    varFrom = "Sequence", 
                                    varTo = "Protein"), 
                 regexp = "No match found between field")
    ## Cannot create a link between an assay and itself
    expect_error(.create_assay_link(rdFrom = rowData(feat1[["psms"]]), 
                                    rdTo = rowData(feat1[["peptides"]]),
                                    from = "psms",
                                    to = "psms",
                                    varFrom = "Sequence", 
                                    varTo = "Sequence"),
                 regexp = "Creating an AssayLink between an assay")
})

test_that(".create_assay_link_multi", {
    ## Artificially duplicate psms assay
    expect_identical(alMulti@name, "peptides")
    expect_identical(alMulti@from, c("psms", "psms2"))
    expect_identical(alMulti@fcol, c("Sequence", "Sequence"))
    ## we use expect_equivalent and not identical  because the `hits` in 
    ## `feat1[["peptides"]]` are sorted by `Sequence` during `aggregateFeatures`
    expect_equivalent(alMulti@hits, List(psms = assayLink(feat1, "peptides")@hits,
                                        psms2 = assayLink(feat1, "peptides")@hits))
    ## Test errors
    expect_error(.create_assay_link_multi(featMultiP, from = c("psms2"), 
                                          to = "peptides",  
                                          varsFrom = c("Sequence", "Sequence"),
                                          varTo = "Sequence"),
                 regexp = "Length of 'from' and length of 'varsFrom'")
})
    
test_that(".update_assay_links", {
    ## Test correct usage
    feat2 <- feat1
    feat2@assayLinks[["peptides"]] <- AssayLink("peptides") # remove AL
    expect_equal(feat1, 
                 .update_assay_links(feat2, feat1@assayLinks[["peptides"]]))
    ## Test errors
    ## `to` assay not present in Feaures object
    alWrong <- AssayLink("foo")
    expect_error(.update_assay_links(feat1, alWrong),
                 regexp = "Assay links names are wrong")
    ## `from`` assay not present in Feaures object
    alWrong <- AssayLink("peptides", "bar")
    expect_error(.update_assay_links(feat1, alWrong),
                 regexp = "@from not valid")
    ## Update AssayLink with an AssayLink object containing wrong rownames
    ## Wrong rownames from the parent assay
    alWrong <- assayLink(feat1, "peptides")
    elementMetadata(alWrong@hits)$names_from <- 
        paste(elementMetadata(alWrong@hits)$names_from, "foo")
    expect_error(.update_assay_links(feat1, alWrong),
                 regexp = "The AssayLink metadata 'names_from'")
    ## Wrong rownames from the child assay 
    alWrong <- assayLink(feat1, "peptides")
    elementMetadata(alWrong@hits)$names_to <- 
        paste(elementMetadata(alWrong@hits)$names_to, "bar")
    expect_error(.update_assay_links(feat1, alWrong),
                 regexp = "The AssayLink metadata 'names_to'")
})

test_that(".update_assay_links_multi_parents", {
    ## Test correct usage
    feat3 <- featMultiP
    feat3@assayLinks@listData$peptides <- alMulti
    expect_identical(feat3, 
                     .update_assay_links_multi_parents(featMultiP, alMulti))
    ## Test errors
    ## Wrong rownames from on of the parent assays
    alMulti@hits$psms@elementMetadata$names_from[1] <- "foo"
    expect_error(.update_assay_links_multi_parents(featMultiP, alMulti), 
                 regexp = "'names_from' does not match")
    ## Wrong rownames from the child assay
    alMulti@hits$psms@elementMetadata$names_from[1] <- "PSM1" ## reset correct
    alMulti@hits$psms@elementMetadata$names_to[1] <- "foo"
    expect_error(.update_assay_links_multi_parents(featMultiP, alMulti), 
                 regexp = "'names_to' does not match")
})


test_that("addAssayLink", {
    feat2 <- feat1
    feat2@assayLinks[["peptides"]] <- AssayLink("peptides") # remove AL
    ## Test correct usage
    expect_equivalent(feat1, 
                      addAssayLink(feat2, from = "psms", to = "peptides",
                                   varFrom = "Sequence", varTo = "Sequence"))
    expect_equivalent(feat1, 
                      addAssayLink(feat2, from = 1, to = 2,
                                   varFrom = "Sequence", varTo = "Sequence"))
    ## Test subsetting still works 
    expect_identical(
        dims(addAssayLink(feat2, from = 1, to = 2, varFrom = "Sequence", 
                          varTo = "Sequence")["PSM1", , ]),
        matrix(1:2, ncol = 1, dimnames = list(NULL, "psms")))
})


test_that("addAssayLinkOneToOne", {
    ## Test correct usage
    psms <- feat1[["psms"]]
    feat3 <- Features(list(psms1 = psms, psms2 = psms), metadata = list())
    feat4 <- addAssayLinkOneToOne(feat3, from = "psms1", to = "psms2")
    al <- assayLink(feat4, "psms2")
    expect_identical(al@name, "psms2")
    expect_identical(al@from, "psms1")
    expect_identical(al@fcol, "OneToOne")
    expect_identical(al@hits@from, 1:10)
    expect_identical(al@hits@to, 1:10)
    # Still works after reordering the rows of one assay
    feat3 <- Features(list(psms1 = psms, psms2 = psms[nrow(psms):1]), 
                      metadata = list())
    feat4 <- addAssayLinkOneToOne(feat3, "psms1", "psms2")
    expect_equivalent(feat3, feat4)
    ## Test subsetting still works 
    expect_identical(feat4["PSM1", , ][["psms1"]],
                     feat4["PSM1", , ][["psms2"]])
    ## Test errors
    ## The assays don't have same size 
    expect_error(addAssayLinkOneToOne(feat1, from = "psms", to = "peptides"),
                 regexp = "must have the same number of rows.")
    ## The rownames don't match
    psms2 <- psms
    rownames(psms2) <- paste(rownames(psms2), "foo")
    feat3 <- Features(list(psms1 = psms, psms2 = psms2), 
                      metadata = list())
    expect_error(addAssayLinkOneToOne(feat3, "psms1", "psms2"),
                 regexp = "Different rownames found")
})

test_that("addAssayLinkMultiParent", {
    ## Test correct use
    expect_identical(alMulti, 
                     addAssayLinkMultiParent(featMultiP, 
                                             from = c("psms", "psms2"), 
                                             to = "peptides", 
                                             varsFrom = c("Sequence", "Sequence"),
                                             varTo = "Sequence")@assayLinks$peptides)
    ## Test warning
    expect_warning(addAssayLinkMultiParent(feat2, 
                                           from = "psms", to = "peptides",
                                           varsFrom = "Sequence", 
                                           varTo = "Sequence"), 
                   regexp = "Only 1 parent supplied, calling 'addAssay'")
    ## Test subsetting still works 
    # expect_identical(
    #     dims(addAssayLinkMultiParent(featMultiP, 
    #                                  from = c("psms", "psms2"), 
    #                                  to = "peptides", 
    #                                  varsFrom = c("Sequence", "Sequence"),
    #                                  varTo = "Sequence")["PSM1", , ]),
    #     matrix(1:2, ncol = 1, dimnames = list(NULL, "psms")))
})
