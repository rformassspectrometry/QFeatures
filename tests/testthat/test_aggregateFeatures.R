data(feat1)
se <- feat1[[1]]

test_that("aggregateFeatures,SummarizedExperiment: errors and message", {
    ## Error: Missing fcol
    expect_error(.aggregateQFeatures(se), regexp = "'fcol' is required")
    ## Error: fcol not present in rowData
    expect_error(.aggregateQFeatures(se, fcol = "missing_fcol"),
                 regexp = "'fcol' not found")
    ## Aggregation with missing data creates message
    data("se_na2")
    ## Message: missing data in assay
    expect_message(.aggregateQFeatures(se_na2, fcol = "nNA", fun = colSums),
                   regexp = "Your quantitative data contain missing values")
    ## Message: missing data in assay AND rowData
    rowData(se_na2)[1, 2] <- NA
    expect_message(.aggregateQFeatures(se_na2, fcol = "nNA", fun = colSums),
                   regexp = "Your quantitative and row data contain missing values")
})

test_that("aggregateFeatures,SummarizedExperiment with 'fun = sum'", {
    aggSE <- .aggregateQFeatures(se, fcol = "Sequence", fun = colSums)

    ## checking quantiation data
    assay1 <- matrix(as.numeric(c(sum(1:3), sum(4:6), sum(7:10),
                                  sum(11:13), sum(14:16), sum(17:20))),
                     ncol = 2,
                     dimnames = list(c("SYGFNAAR", "ELGNDAYK", "IAEESNFPFIK"),
                                     c("S1", "S2")))
    assay1 <- assay1[levels(factor(rowData(feat1[[1]])$Sequence)), ]
    expect_identical(assay1, assay(aggSE))

    ## checking rowData
    Sequence <- rownames(assay1)
    Protein <- c("ProtA", "ProtB", "ProtA")
    .n <- c(3L, 4L, 3L)
    location <- c("Mitochondrion", "unknown", "Mitochondrion")
    coldat1 <- DataFrame(Sequence = Sequence,
                         Protein = Protein,
                         location = location,
                         .n = .n,
                         row.names = rownames(assay1))
    expect_identical(coldat1, rowData(aggSE))
})

test_that("aggregateFeatures,SummarizedExperiment with 'fun = median'", {
    aggSE <- .aggregateQFeatures(se, fcol = "Sequence",
                                fun = matrixStats::colMedians)

    ## checking quantiation data
    assay1 <- matrix(as.numeric(c(median(1:3), median(4:6), median(7:10),
                                  median(11:13), median(14:16), median(17:20))),
                     ncol = 2,
                     dimnames = list(c("SYGFNAAR", "ELGNDAYK", "IAEESNFPFIK"),
                                     c("S1", "S2")))
    assay1 <- assay1[order(rownames(assay1)), ]
    expect_identical(assay1, assay(aggSE))

    ## checking rowData
    Sequence <- rownames(assay1)
    Protein <- c("ProtA", "ProtB", "ProtA")
    .n <- c(3L, 4L, 3L)
    location <- c("Mitochondrion", "unknown", "Mitochondrion")
    coldat1 <- DataFrame(Sequence = Sequence,
                         Protein = Protein,
                         location = location,
                         .n = .n,
                         row.names = rownames(assay1))
    expect_equal(coldat1, rowData(aggSE))
})

test_that("aggregateFeatures,SummarizedExperiment return class (issue 78)", {
    library(SingleCellExperiment)
    ## Aggregating an SEs produces an SE
    expect_identical(class(aggregateFeatures(se, fcol = "Sequence"))[[1]],
                     "SummarizedExperiment")
    ## Aggregating an SCEs produces an SCE
    expect_identical(class(aggregateFeatures(as(se, "SingleCellExperiment"),
                                             fcol = "Sequence"))[[1]],
                     "SingleCellExperiment")
})

test_that("aggregateFeatures,QFeatures: empty and errors", {
    expect_identical(QFeatures(), aggregateFeatures(QFeatures()))
    expect_error(aggregateFeatures(feat1, name = "psms"),
                 regexp = "There's already an assay named 'psms'")
})

test_that("aggregateFeatures,QFeatures: check links and subsetting", {
    feat1 <- aggregateFeatures(feat1, "psms", fcol = "Sequence",
                               name = "peptides", fun = colSums)
    ## Checking the aggregated assay is correctly added
    expect_identical(names(feat1), c("psms", "peptides"))
    expect_identical(dims(feat1),
                     matrix(c(10L, 2L, 3L, 2L), ncol = 2,
                            dimnames = list(NULL, c("psms", "peptides"))))

    ## Checking assayLinks
    alink <- feat1@assayLinks[[2]]
    expect_identical(alink@from, "psms")
    expect_identical(alink@fcol, "Sequence")
    hits1 <- Hits(from = 1:10,
                  to = c(rep(3, 3), rep(1, 3), rep(2, 4)),
                  names_from = paste0("PSM", 1:10),
                  names_to = c(rep("SYGFNAAR", 3),
                               rep("ELGNDAYK", 3),
                               rep("IAEESNFPFIK", 4)),
                  nLnode = 10L,
                  nRnode = 3L,
                  sort.by.query = TRUE)
    expect_identical(alink@hits, hits1)

    ## Checking subsetting still works
    featsub <- feat1["IAEESNFPFIK", , ]
    expect_identical(dims(featsub),
                     matrix(c(4L, 2L, 1L, 2L), ncol = 2,
                            dimnames = list(NULL, c("psms", "peptides"))))
    ## The rowData should contain only the Sequence used for subsetting
    expect_identical(unique(rowData(featsub[["peptides"]])$Sequence),
                     "IAEESNFPFIK")
})

test_that("aggregateFeatures,QFeatures: aggcounts", {
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
