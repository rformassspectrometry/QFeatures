data(feat1)
se <- feat1[[1]]

test_that("aggregateFeatures,SummarizedExperiment: errors and message", {
    ## Error: Missing fcol
    expect_error(.aggregateQFeatures(se), regexp = "'fcol' is required")
    ## Error: fcol not present in rowData
    expect_error(.aggregateQFeatures(se, fcol = "missing_fcol"),
                 regexp = "'fcol' not found")
    ## Error: fcol is a list
    rowData(se)$lseq <- as.list(rowData(se)$Sequence)
    expect_error(.aggregateQFeatures(se, fcol = "lseq", fun = colSums),
                 "'fcol' must refer to an atomic vector or a sparse matrix")

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
                 regexp = "one or more assays named: 'psms'")
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


test_that("aggregate by matrix and (atomic) vector work (1)", {
    ## only features that map uniquely, aggregation should thus be
    ## identical whether we use a vector or an adjacency matrix
    se <- feat1[[1]]
    ## ----------------------------------------------------
    ## (1) PSM to peptide aggregation
    ## generated with PSMatch::makeAdjacencyMatrix(rowData(se)$Sequence)
    adjSequence <- structure(
        c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
          1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
        .Dim = c(10L, 3L),
        .Dimnames = list(
            paste0("PSM", 1:10),
            c("SYGFNAAR", "ELGNDAYK", "IAEESNFPFIK")))
    ## 1.1 aggregate PSMs to peptides by vector
    se1 <- aggregateFeatures(se, "Sequence", colSums)
    ## 1.2 aggregate PSMs to peptides by matrix
    rowData(se)$adjacencyMatrix <- adjSequence
    expect_error(aggregateFeatures(se, "adjacencyMatrix",
                                   MsCoreUtils::colSumsMat),
                 "'fcol' must refer to an atomic vector or a sparse matrix")
    ## Error when adjancencyMatrix already present
    expect_error(adjacencyMatrix(se) <- adjSequence,
                 "Found an existing variable adjacencyMatrix.")
    rowData(se)[["adjacencyMatrix"]] <- NULL
    adjacencyMatrix(se) <- adjSequence
    se2 <- aggregateFeatures(se, "adjacencyMatrix", MsCoreUtils::colSumsMat)
    ## order of rows isn't necessary the same
    rnms <- rownames(se1)
    expect_identical(assay(se1)[rnms, ], assay(se2)[rnms, ])
    expect_identical(colData(se1), colData(se2))
    ## below not identical/equal because '.n' is named in se2
    expect_equivalent(rowData(se1)[rnms, ], rowData(se2)[rnms, ])
    ## ----------------------------------------------------
    ## (2) Peptide to protein aggregation
    ## generated with PSMatch::makeAdjacencyMatrix(rowData(se2)$Protein)
    adjProtein <- structure(
        c(1, 1, 0, 0, 0, 1),
        .Dim = 3:2,
        .Dimnames = list(c("SYGFNAAR", "ELGNDAYK", "IAEESNFPFIK"),
                         c("ProtA", "ProtB")))
    ## 2.1 aggregate peptides to proteins by (atomic) vector
    se3 <- aggregateFeatures(se1, "Protein", colSums)
    ## 2.2 aggregate peptides to proteins by matrix
    adjacencyMatrix(se2) <- adjProtein
    se4 <- aggregateFeatures(se2, "adjacencyMatrix", MsCoreUtils::colSumsMat)
    ## order of rows isn't necessary the same
    rnms <- rownames(se3)
    expect_identical(assay(se3)[rnms, ], assay(se4)[rnms, ])
    expect_identical(colData(se3), colData(se4))
    ## below not identical/equal because '.n' is named in se2
    expect_equivalent(rowData(se3)[rnms, ], rowData(se4)[rnms, ])
})


test_that("aggregate by matrix and (atomic) vector work (2)", {
    ## Change last PSM/peptide to be shared among proteins B and C
    se <- feat1[[1]]
    rowData(se)[10, "Sequence"] <- "PEPTIDE"
    rowData(se)[10, "Protein"] <- "ProtB;ProtC"
    ## prepare data
    se <- aggregateFeatures(se, "Sequence", colSums)
    ## -----------------------------------------------------
    ## aggregate peptides to proteins by vector
    se1 <- aggregateFeatures(se, "Protein", colSums)
    ## generated with PSMatch::makeAdjacencyMatrix(rowData(se)$Protein)
    adjProtein <- structure(
        c(1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0),
        .Dim = 4:3,
        .Dimnames = list(c("ELGNDAYK", "IAEESNFPFIK", "PEPTIDE", "SYGFNAAR"),
                         c("ProtA", "ProtB", "ProtC")))
    adjacencyMatrix(se) <- adjProtein
    se2 <- aggregateFeatures(se, "adjacencyMatrix", MsCoreUtils::colSumsMat)
    ## only ProtA is composed of unique peptides only - only protein
    ## with identical aggregation results
    expect_identical(assay(se1)["ProtA", ], assay(se2)["ProtA", ])
    expect_identical(colData(se1), colData(se2))
    k <- intersect(names(rowData(se1)), names(rowData(se2)))
    ## below not identical/equal because '.n' is named in se2
    expect_equivalent(rowData(se1)["ProtA", k], rowData(se2)["ProtA", k])
})

test_that("aggregateFeatures,QFeatures: aggregate multiple assays", {
    data("feat3")
    expect_warning(feat3 <- feat3[, , 1:3],
                   regexp = "experiments' dropped; see 'drops")
    ii <- names(feat3)
    feat3aggr <- aggregateFeatures(feat3, i = ii,
                                   fcol = rep("Protein", 3),
                                   name = paste0("prots", 1:3),
                                   fun = colSums)
    ## Alternatively there is no need to repeat fcol 3x
    expect_identical(feat3aggr,
                     aggregateFeatures(feat3, i = ii, fcol = "Protein",
                                       name = paste0("prots", 1:3),
                                       fun = colSums))
    ## Checking the aggregated assay is correctly added
    expect_identical(dims(feat3aggr),
                     matrix(c(7L, 8L, 10L, rep(2L, 3), rep(c(2L, 2L, 4L), 2)),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(NULL, c(ii, paste0("prots", 1:3)))))
    ## Checking subsetting still works
    featsub <- feat3aggr["ProtA", , ]
    expect_identical(dims(featsub),
                     matrix(c(6L, 4L, 6L, rep(1L, 3), rep(c(2L, 2L, 4L), 2)),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(NULL, c(ii, paste0("prots", 1:3)))))
    ## The rowData should contain only the Sequence used for subsetting
    expect_identical(unique(unlist(rownames(featsub)[4:6])),
                     "ProtA")
    ## Test errors
    ## One assay is already present
    expect_error(aggregateFeatures(feat3, i = ii, fcol = rep("Protein", 3),
                                   name = paste0("psms", 1:3),
                                   fun = colSums),
                 regexp = "named: 'psms1', 'psms2'")
    ## 'i' and 'name' must have same length
    expect_error(aggregateFeatures(feat3, i = ii, fcol = rep("Protein", 3),
                                   name = "prot1",
                                   fun = colSums),
                 regexp = "'i' and 'name' must have same length")
    ## 'i' and 'fcol' must have same length
    expect_error(aggregateFeatures(feat3, i = ii, fcol = rep("Protein", 4),
                                   name = paste0("prot", 1:3),
                                   fun = colSums),
                 regexp = "'i' and 'fcol' must have same length")
})

test_that("aggregateFeatures,QFeatures: invariant columns discarded", {
    data("feat3")
    expect_warning(feat3 <- feat3[, , 1:3],
                   regexp = "experiments' dropped; see 'drops")
    ii <- names(feat3)
    for (i in ii) {
       rowData(feat3[[i]])$variant <- seq_along(nrow(rowData(feat3[[i]])))
    }
    rowData(feat3[[1]])$psms1col <- 1

   feat3aggr <- aggregateFeatures(feat3, i = ii, fcol = "Protein",
                                       name = paste0("prots", 1:3),
                                       fun = colSums)
    ## Checking the aggregated assay rowData doesn't contain the variant columns
    expect_identical(dims(rowData(feat3aggr)),
                     matrix(c(7L, 7L, 8L, 6L, 10L, 6L, rep(c(2L, 4L), 3)),
                            nrow = 6, byrow = TRUE,
                            dimnames = list(c(ii, paste0("prots", 1:3)), NULL)))
})
