data(feat1)
feat1 <- aggregateFeatures(feat1, "psms", "Sequence", "peptides")
## Artificially duplicate psms assay to test multi-parent linking
featMultiP <- addAssay(feat1, feat1[["psms"]], name = "psms2")
alMulti <- .create_assay_link(featMultiP, 
                              from = c("psms", "psms2"),
                              to = "peptides",  
                              varFrom = c("Sequence", "Sequence"),
                              varTo = "Sequence")

test_that("empty AssayLink", {
    ## Test the constructor for generating an empty AssayLink object
    al <- AssayLink("foo")
    naChar <- character(1)
    naChar[1] <- NA ## Get an NA of class "character
    expect_true(validObject(al))
    expect_identical(al@name, "foo")
    expect_identical(al@from, naChar)
    expect_identical(al@hits, Hits())
    expect_null(show(al))
})

test_that("empty AssayLinks", {
    ## Test the constructor for generating an empty AssayLinks object
    als <- AssayLinks()
    expect_true(validObject(als))
    expect_identical(length(als), 0L)
    expect_null(show(als))
})

test_that("Subset AssayLinks", {
    als <- assayLinks(feat1, "peptides")
    als <- als[list(peptides = "ELGNDAYK", psms = paste0("PSM", 4:6))]
    ## Check that only the desired peptide remains
    expect_identical("ELGNDAYK", unique(elementMetadata(als[["peptides"]]@hits)$names_to))
    ## Check that only the desired PSMS remains
    expect_identical(paste0("PSM", 4:6), unique(elementMetadata(als[["peptides"]]@hits)$names_from))
})

test_that(".getHits", {
    ## Test one to one matching, same ordering
    hits <- .get_Hits(rdFrom = DataFrame(A = 1:10, row.names = letters[1:10]),
                      rdTo = DataFrame(B = 1:10, row.names = LETTERS[1:10]),
                      varFrom = "A", varTo = "B")
    ## Test the expected number of hits
    expect_identical(length(hits), 10L) 
    ## Test the rownames are correctly matched
    expect_identical(elementMetadata(hits)$names_from, letters[1:10]) 
    expect_identical(elementMetadata(hits)$names_to, LETTERS[1:10]) 
    
    ## Test one to one matching, different ordering
    hits <- .get_Hits(rdFrom = DataFrame(A = 1:10, row.names = letters[1:10]),
                      rdTo = DataFrame(B = 10:1, row.names = LETTERS[10:1]),
                      varFrom = "A", varTo = "B")
    ## Test the expected number of hits
    expect_identical(length(hits), 10L) 
    ## Test the rownames are correctly matched
    expect_identical(elementMetadata(hits)$names_from, letters[1:10]) 
    expect_identical(elementMetadata(hits)$names_to, LETTERS[1:10]) 
    
    ## Test many to one matching 
    hits <- .get_Hits(rdFrom = DataFrame(A = rep(1:10, each = 2), row.names = letters[1:20]),
                      rdTo = DataFrame(B = 1:10, row.names = LETTERS[1:10]),
                      varFrom = "A", varTo = "B")
    ## Test the expected number of hits
    expect_identical(length(hits), 20L) 
    ## Test the rownames are correctly matched
    expect_identical(elementMetadata(hits)$names_from, letters[1:20]) 
    expect_identical(elementMetadata(hits)$names_to, rep(LETTERS[1:10], each = 2)) 
    
    ## Test many to one matching, with some features missing
    hits <- .get_Hits(rdFrom = DataFrame(A = rep(1:10, each = 2), row.names = letters[1:20]),
                      rdTo = DataFrame(B = 1:5, row.names = LETTERS[1:5]),
                      varFrom = "A", varTo = "B")
    ## Test the expected number of hits
    expect_identical(length(hits), 10L) 
    ## Test the rownames are correctly matched
    expect_identical(elementMetadata(hits)$names_from, letters[1:10]) 
    expect_identical(elementMetadata(hits)$names_to, rep(LETTERS[1:5], each = 2)) 
    
    ## Test errors 
    ## No match found 
    expect_error(.get_Hits(rdFrom = DataFrame(A = 1:5),  rdTo = DataFrame(B = 6:10),
                           varFrom = "A", varTo = "B"), 
                 regexp = "No match found")
})

test_that(".create_assay_link", {
    ## Test 1 parent to 1 child link
    ## we use expect_equivalent and not identical  because the `hits` in 
    ## `feat1[["peptides"]]` are sorted by `Sequence` during `aggregateFeatures`
    expect_equivalent(feat1@assayLinks[[2]], 
                      .create_assay_link(feat1, from = "psms", to = "peptides",
                                         varFrom = "Sequence", varTo = "Sequence"))
    ## Test one to one link based on rownaes
    al <- .create_assay_link(featMultiP, from = "psms", to = "psms2")
    expect_identical(al@from, "psms")
    expect_identical(al@fcol, "._rownames")
    expect_identical(elementMetadata(al@hits)$names_from, paste0("PSM", 1:10))
    expect_identical(elementMetadata(al@hits)$names_to, paste0("PSM", 1:10))
    ## Test 2 parents to 1 child
    expect_identical(alMulti@name, "peptides")
    expect_identical(alMulti@from, c("psms", "psms2"))
    expect_identical(alMulti@fcol, c("Sequence", "Sequence"))
    expect_true(inherits(alMulti@hits, "List"))
    expect_identical(alMulti@hits[[1]], alMulti@hits[[2]])
    expect_equivalent(alMulti@hits[[1]], feat1@assayLinks[["peptides"]]@hits)
    
    ## Test errors
    ## Cannot create a link between an assay and itself
    expect_error(.create_assay_link(feat1, from = "psms", to = "psms",
                                    varFrom = "Sequence", varTo = "Sequence"),
                 regexp = "itself is not allowed")
    ## The number of `from` assays and the number of `varFrom` do not match
    expect_error(.create_assay_link(feat1, from = "psms", to = "peptides",  
                                    varFrom = c("Sequence", "Sequence"),
                                    varTo = "Sequence"),
                 regexp = "must have same length")
})

    
test_that(".update_assay_links", {
    ## Test 1 parent to 1 child link
    feat2 <- feat1
    feat2@assayLinks[["peptides"]] <- AssayLink("peptides") # remove AL
    expect_equal(feat1, 
                 .update_assay_links(feat2, feat1@assayLinks[["peptides"]]))
    ## Test 2 parents to 1 child link
    feat3 <- .update_assay_links(featMultiP, alMulti)
    expect_identical(assayLink(feat3,"peptides")@from, alMulti@from)
    expect_identical(assayLink(feat3,"peptides")@fcol, alMulti@fcol)
    expect_identical(assayLink(feat3,"peptides")@hits, alMulti@hits)
    
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
                 regexp = "does not match the rownames")
    ## Same but with multiple-parent link
    ## Wrong rownames from on of the parent assays
    alMulti@hits$psms@elementMetadata$names_from[1] <- "foo"
    expect_error(.update_assay_links(featMultiP, alMulti), 
                 regexp = "'names_from' does not match")
    ## Wrong rownames from the child assay
    alMulti@hits$psms@elementMetadata$names_from[1] <- "PSM1" ## reset correct
    alMulti@hits$psms@elementMetadata$names_to[1] <- "foo"
    expect_error(.update_assay_links(featMultiP, alMulti), 
                 regexp = "'names_to' does not match")
})

test_that("addAssayLink", {
    feat2 <- feat1
    feat2@assayLinks[["peptides"]] <- AssayLink("peptides") # remove AL
    ## Test 1 parent to 1 child link
    expect_identical(addAssayLink(feat2, from = "psms", to = "peptides",
                                  varFrom = "Sequence", varTo = "Sequence"),
                     addAssayLink(feat2, from = 1, to = 2,
                                  varFrom = "Sequence", varTo = "Sequence"))
    expect_equivalent(feat1, 
                      addAssayLink(feat2, from = "psms", to = "peptides",
                                   varFrom = "Sequence", varTo = "Sequence"))
    ## Test 2 parents to 1 child link
    expect_identical(alMulti, 
                     addAssayLink(featMultiP, 
                                  from = c("psms", "psms2"), 
                                  to = "peptides", 
                                  varFrom = c("Sequence", "Sequence"),
                                  varTo = "Sequence")@assayLinks$peptides)
    expect_identical(alMulti, 
                     addAssayLink(featMultiP, 
                                  from = c(1, 3), 
                                  to = "peptides", 
                                  varFrom = c("Sequence", "Sequence"),
                                  varTo = "Sequence")@assayLinks$peptides)
})


test_that("addAssayLinkOneToOne", {
    ## Test correct usage
    psms <- feat1[["psms"]]
    feat2 <- Features(list(psms1 = psms, psms2 = psms))
    feat2 <- addAssayLinkOneToOne(feat2, from = "psms1", to = "psms2")
    al <- assayLink(feat2, "psms2")
    expect_identical(al@name, "psms2")
    expect_identical(al@from, "psms1")
    expect_identical(al@fcol, "._oneToOne")
    expect_identical(al@hits@from, 1:10)
    expect_identical(al@hits@to, 1:10)
    # Still works after reordering the rows of one assay
    feat3 <- Features(list(psms1 = psms, psms2 = psms[nrow(psms):1]))
    feat3 <- addAssayLinkOneToOne(feat3, "psms1", "psms2")
    expect_equivalent(feat2@assayLinks[[2]], feat3@assayLinks[[2]])
    
    ## Test errors
    ## The assays don't have same size 
    expect_error(addAssayLinkOneToOne(feat1, from = "psms", to = "peptides"),
                 regexp = "must have the same number of rows.")
    ## The rownames don't match
    psms2 <- psms
    rownames(psms2) <- paste(rownames(psms2), "foo")
    feat2 <- Features(list(psms1 = psms, psms2 = psms2))
    expect_error(addAssayLinkOneToOne(feat2, "psms1", "psms2"),
                 regexp = "Different rownames found")
    ## One to one link is not allowed between an assay and itself
    expect_error(addAssayLinkOneToOne(feat2, "psms1", "psms1"),
                 regexp = "itself is not allowed")
    ## Adding one to one links is not support for multiple parents 
    expect_error(addAssayLinkOneToOne(feat2, c("psms1", "psms2"), "psms3"),
                 regexp = "One to one links are not supported for multiple parents.")
})


