data(feat1)
psms <- feat1[["psms"]]
feat1 <- aggregateFeatures(feat1, "psms", "Sequence", "peptides")
feat2 <- addAssay(feat1[, , 1], feat1[["peptides"]], name = "peptides")

test_that(".createAssayLink", {
    ## Correct use
    ## we here use expect_equivalent and not identical  because the `hits` in 
    ## `feat1` are sorted by `Sequence` by `aggregateFeatures`
    expect_equivalent(feat1@assayLinks[[2]], 
                      .createAssayLink(seFrom = feat2[["psms"]], 
                                       seTo = feat2[["peptides"]],
                                       nameFrom = "psms",
                                       nameTo = "peptide",
                                       varFrom = "Sequence", 
                                       varTo = "Sequence"))
    ## No match between feature variables
    expect_error(.createAssayLink(seFrom = feat2[["psms"]], 
                                  seTo = feat2[["peptides"]],
                                  nameFrom = "psms", nameTo = "peptide",
                                  varFrom = "Sequence", varTo = "Protein"), 
                 regexp = "No match found between field")
    ## Cannot create a link between an assay and itself
    expect_error(.createAssayLink(seFrom = feat2[["psms"]], 
                                  seTo = feat2[["psms"]],
                                  nameFrom = "A", nameTo = "A",
                                  varFrom = "Sequence", varTo = "Sequence"),
                 regexp = "Creating an AssayLink between an assay")
})

test_that(".createAssayLinkOneToOne", {
    al <- .createAssayLinkOneToOne(seFrom = psms, seTo = psms,
                                   nameFrom = "A", nameTo = "B", 
                                   varComm = "Var")
    expect_identical(al@name, "B")
    expect_identical(al@from, "A")
    expect_identical(al@fcol, "Var")
    expect_identical(al@hits@from, 1:10)
    expect_identical(al@hits@to, 1:10)
    al@fcol <- "oneToOneID"
    expect_identical(al, .createAssayLinkOneToOne(seFrom = psms, seTo = psms,
                                                  nameFrom = "A", nameTo = "B"))
    # Not the same dimensions 
    expect_error(.createAssayLinkOneToOne(seFrom = psms, 
                                          seTo = feat2[["peptides"]],
                                          nameFrom = "A", nameTo = "B", 
                                          varComm = "Var"),
                 regexp = "must have the same number of rows")
})

test_that(".addAssayLink", {
    feat2noAL <- feat2
    feat2noAL@assayLinks[["peptides"]] <- AssayLink("peptides")
    expect_equal(feat2, 
                 .addAssayLink(feat2noAL, feat2@assayLinks[["peptides"]]))
})


test_that("createAssayLink", {
    expect_equivalent(feat1, 
                      createAssayLink(feat2, from = "psms", to = "peptides",
                                      varFrom = "Sequence", varTo = "Sequence"))
    expect_equivalent(feat1, 
                      createAssayLink(feat2, from = 1, to = 2,
                                      varFrom = "Sequence", varTo = "Sequence"))
    expect_equivalent(feat1, 
                      createAssayLink(feat2, from = 1, to = 2, 
                                      varFrom = "Sequence"))
})


test_that("createAssayLinkOneToOne", {
    feat3 <- Features(list(psms1 = psms, psms2 = psms), metadata = list())
    feat4 <- createAssayLinkOneToOne(feat3, from = "psms1", to = "psms2", 
                                     varCommon = "Var")
    al <- assayLinks(feat4, "psms2")[[1]]
    expect_identical(al@name, "psms2")
    expect_identical(al@from, "psms1")
    expect_identical(al@fcol, "Var")
    expect_identical(al@hits@from, 1:10)
    expect_identical(al@hits@to, 1:10)
    feat4@assayLinks[["psms2"]]@fcol <- "oneToOneID"
    expect_identical(feat4, 
                     createAssayLinkOneToOne(feat3, from = "psms1", to = "psms2"))
})
