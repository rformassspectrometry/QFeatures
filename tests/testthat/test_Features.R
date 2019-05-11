data(feat1)

test_that("empty Features", {
    feat0 <- Features()
    expect_true(validObject(feat0))
    expect_true(isEmpty(feat0))
})


test_that("Manual Features", {
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
    feat2 <- Features(list(psms = psms), colData = coldata)
    expect_true(validObject(feat2))
    ## load data
    ## data(feat1)
    ## expect_equivalent(feat1, feat2)
})


test_that("combineFeatures(fun = sum)", {
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

    hits1 <- Hits(from = c(rep(1, 3), rep(2, 4), rep(3, 3)),
                  to = c(4:10, 1:3),
                  names_to = c(rep("ELGNDAYK", 3),
                               rep("IAEESNFPFIK", 4),
                               rep("SYGFNAAR", 3)),                      
                  nLnode = 3L,
                  nRnode = 10L,
                  names_from = paste0("PSM", c(4:10, 1:3)),
                  sort.by.query = TRUE)
    expect_identical(alink@hits, hits1)    
})

test_that("combineFeatures(fun = median)", {
    feat1 <- combineFeatures(feat1, "psms", fcol = "Sequence", name = "peptides", fun = median)
    expect_identical(names(feat1), c("psms", "peptides"))

    ## checking quantiation data
    assay1 <- matrix(c(median(1:3), median(4:6), median(7:10),
                       median(11:13), median(14:16), median(17:20)),
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

    hits1 <- Hits(from = c(rep(1, 3), rep(2, 4), rep(3, 3)),
                  to = c(4:10, 1:3),
                  names_to = c(rep("ELGNDAYK", 3),
                               rep("IAEESNFPFIK", 4),
                               rep("SYGFNAAR", 3)),                      
                  nLnode = 3L,
                  nRnode = 10L,
                  names_from = paste0("PSM", c(4:10, 1:3)),
                  sort.by.query = TRUE)
    expect_identical(alink@hits, hits1)    
})

test_that("addAssay", {
    data(feat1)
    assay2 <- feat1[[1]]
    feat1 <- addAssay(feat1, assay2, name = "psms2")
    expect_identical(names(feat1), c("psms", "psms2"))
    expect_identical(feat1[[1]], feat1[[2]])
})
