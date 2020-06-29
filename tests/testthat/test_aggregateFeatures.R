data(feat1)

test_that("aggregateFeatures: empty and errors", {
    expect_identical(Features(), aggregateFeatures(Features()))
    expect_error(aggregateFeatures(feat1, name = "psms"))
    expect_error(aggregateFeatures(feat1))
    expect_error(aggregateFeatures(feat1, fcol = "missing_fcol"))
})

test_that("aggregateFeatures(fun = sum)", {
    feat1 <- aggregateFeatures(feat1, "psms", fcol = "Sequence",
                               name = "peptides", fun = colSums)
    expect_identical(names(feat1), c("psms", "peptides"))

    ## checking quantiation data
    assay1 <- matrix(c(sum(1:3), sum(4:6), sum(7:10),
                       sum(11:13), sum(14:16), sum(17:20)),
                     ncol = 2,
                     dimnames = list(c("SYGFNAAR", "ELGNDAYK", "IAEESNFPFIK"),
                                     c("S1", "S2")))
    assay1 <- assay1[levels(factor(rowData(feat1[[1]])$Sequence)), ]
    expect_equal(assay1, assay(feat1, "peptides"))

    ## checking rowData
    Sequence <- rownames(assay1)
    Protein <- c("ProtA", "ProtB", "ProtA")
    .n <- c(3, 4, 3)
    location <- c("Mitochondrion", "unknown", "Mitochondrion")

    coldat1 <- DataFrame(Sequence = Sequence,
                         Protein = Protein,
                         location = location,
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

test_that("aggregateFeatures(fun = median)", {
    data(feat1)
    feat1 <- aggregateFeatures(feat1, "psms", fcol = "Sequence",
                               name = "peptides",
                               fun = matrixStats::colMedians)
    expect_identical(names(feat1), c("psms", "peptides"))

    expect_equal(dims(feat1),
                 matrix(c(10, 2, 3, 2), ncol = 2,
                        dimnames = list(NULL,
                                        c("psms", "peptides"))))

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
    location <- c("Mitochondrion", "unknown", "Mitochondrion")

    coldat1 <- DataFrame(Sequence = Sequence,
                         Protein = Protein,
                         location = location,
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

test_that("aggregateFeatures: aggcounts", {
    data(ft_na)
    ## missing values are propagated
    ft2 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
    se <- ft2[[2]]
    ## assay
    m <- matrix(c(NA, NA, NA, 14, 20, 22), nrow = 2)
    rownames(m) <- 1:2
    colnames(m) <- LETTERS[1:3]
    expect_identical(assay(se), m)
    ## aggcounts
    m <- matrix(c(1, 1, 1, 2, 2, 2), nrow = 2)
    rownames(m) <- 1:2
    colnames(m) <- LETTERS[1:3]
    expect_identical(aggcounts(se), m)
    ## .n variable
    expect_identical(rowData(se)$.n , c(2L, 2L))
})
