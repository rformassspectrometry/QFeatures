data(feat1)

test_that("combineFeatures", {
    feat1 <- combineFeatures(feat1, "psms", fcol = "Sequence", name = "peptides", fun = sum)
    expect_identical(names(feat1), c("psms", "peptides"))

    ## checking quantiation data
    assay1 <- matrix(c(sum(1:3), sum(4:6), sum(7:10),
                       sum(11:13), sum(14:16), sum(17:20)),
                     ncol = 2,
                     dimnames = list(c("SYGFNAAR", "ELGNDAYK", "IAEESNFPFIK"),
                                     c("S1", "S2")))
    assay1 <- assay1[order(rownames(assay1)), ]
    expect_identical(assay1, assay(feat1, "peptides"))

    ## checking rowData
    Sequence <- rownames(assay1)
    Protein <- c("ProtA", "ProtB", "ProtA")
    .n <- c(3, 4, 3)
    names(.n) <- names(Protein) <- names(Sequence) <- Sequence

    coldat1 <- DataFrame(Sequence = Sequence,
                         Protein = Protein,
                         .n = .n,                         
                         row.names = rownames(assay1))
    expect_equal(coldat1, rowData(feat1[["peptides"]]))


    ## checking assayLinks
    alink <- feat1@assayLinks[[2]]
    expect_identical(alink@from, "psms")
    expect_identical(alink@fcol, "Sequence")
    
})
