data(feat2)

se1 <- feat2[[1]]
se2 <- feat2[[2]]
se3 <- feat2[[3]]


test_that("merge SEs (1,2)", {
    se <- QFeatures:::mergeSElist(list(se1, se2))
    expect_identical(ncol(se), ncol(se1) + ncol(se2))
    expect_identical(nrow(se), length(union(rownames(se1), rownames(se2))))
    expect_identical(names(rowData(se)), "Prot")
    ## test assay data
    na2 <- c(S5 = NA, S6 = NA, S7 = NA, S8 = NA)
    na1 <- c(S1 = NA, S2 = NA, S3 = NA, S4 = NA)
    expect_identical(assay(se)["a", ], c(assay(se1)["a", ], na2))
    expect_identical(assay(se)["b", ], c(assay(se1)["b", ], na2))
    expect_identical(assay(se)["c", ], c(assay(se1)["c", ], na2))
    expect_identical(assay(se)["d", ], c(assay(se1)["d", ], na2))
    expect_identical(assay(se)["e", ], c(assay(se1)["e", ], na2))
    expect_identical(assay(se)["f", ], c(assay(se1)["f", ], na2))
    expect_identical(assay(se)["g", ], c(assay(se1)["g", ], na2))
    expect_identical(assay(se)["h", ], c(assay(se1)["h", ], assay(se2)["h", ]))
    expect_identical(assay(se)["i", ], c(assay(se1)["i", ], assay(se2)["i", ]))
    expect_identical(assay(se)["k", ], c(na1, assay(se2)["k", ]))
    ## test rowData
    expect_identical(rowData(se)$Prot, paste0("P", rownames(se)))
    expect_identical(sort(rowData(se)$Prot),
                     sort(union(rowData(se1)$Prot, rowData(se2)$Prot)))
})

test_that("merge SEs (2, 3)", {
    se <- QFeatures:::mergeSElist(list(se2, se3))
    expect_identical(ncol(se), ncol(se2) + ncol(se3))
    expect_identical(nrow(se), length(union(rownames(se2), rownames(se3))))
    expect_identical(names(rowData(se)), "Prot")
    ## test assay data
    na2 <- c(S5 = NA, S6 = NA, S7 = NA, S8 = NA)
    na3 <- c(S9 = NA, S10 = NA, S11 = NA, S12 = NA)
    expect_identical(assay(se)["a", ], c(na2, assay(se3)["a", ]))
    expect_identical(assay(se)["b", ], c(na2, assay(se3)["b", ]))
    expect_identical(assay(se)["l", ], c(na2, assay(se3)["l", ]))
    expect_identical(assay(se)["m", ], c(na2, assay(se3)["m", ]))
    expect_identical(assay(se)["n", ], c(na2, assay(se3)["n", ]))
    expect_identical(assay(se)["h", ], c(assay(se2)["h", ], na3))
    expect_identical(assay(se)["i", ], c(assay(se2)["i", ], na3))
    expect_identical(assay(se)["j", ], c(assay(se2)["j", ], assay(se3)["j", ]))
    expect_identical(assay(se)["k", ], c(assay(se2)["k", ], assay(se3)["k", ]))
    ## test rowData
    expect_identical(rowData(se)$Prot, paste0("P", rownames(se)))
    expect_identical(sort(rowData(se)$Prot),
                     sort(union(rowData(se2)$Prot, rowData(se3)$Prot)))
})

test_that("merge SEs (3, 2)", {
    se <- QFeatures:::mergeSElist(list(se3, se2))
    expect_identical(ncol(se), ncol(se2) + ncol(se3))
    expect_identical(nrow(se), length(union(rownames(se2), rownames(se3))))
    expect_identical(names(rowData(se)), "Prot")
    ## test assay data
    na2 <- c(S5 = NA, S6 = NA, S7 = NA, S8 = NA)
    na3 <- c(S9 = NA, S10 = NA, S11 = NA, S12 = NA)
    expect_identical(assay(se)["a", ], c(assay(se3)["a", ], na2))
    expect_identical(assay(se)["b", ], c(assay(se3)["b", ], na2))
    expect_identical(assay(se)["l", ], c(assay(se3)["l", ], na2))
    expect_identical(assay(se)["m", ], c(assay(se3)["m", ], na2))
    expect_identical(assay(se)["n", ], c(assay(se3)["n", ], na2))
    expect_identical(assay(se)["h", ], c(na3, assay(se2)["h", ]))
    expect_identical(assay(se)["i", ], c(na3, assay(se2)["i", ]))
    expect_identical(assay(se)["j", ], c(assay(se3)["j", ], assay(se2)["j", ]))
    expect_identical(assay(se)["k", ], c(assay(se3)["k", ], assay(se2)["k", ]))
    ## test rowData
    expect_identical(rowData(se)$Prot, paste0("P", rownames(se)))
    expect_identical(sort(rowData(se)$Prot),
                     sort(union(rowData(se2)$Prot, rowData(se3)$Prot)))
})


test_that("merge SEs (1, 2, 3)", {
    se <- QFeatures:::mergeSElist(list(se1, se2, se3))
    expect_identical(ncol(se), ncol(se1) + ncol(se2) + ncol(se3))
    expect_identical(nrow(se),
                     length(Reduce(union, list(rownames(se1), rownames(se2), rownames(se3)))))
    expect_identical(names(rowData(se)), "Prot")
    ## test assay data
    na1 <- c(S1 = NA, S2 = NA, S3 = NA, S4 = NA)
    na2 <- c(S5 = NA, S6 = NA, S7 = NA, S8 = NA)
    na3 <- c(S9 = NA, S10 = NA, S11 = NA, S12 = NA)
    expect_identical(assay(se)["a", ], c(assay(se1)["a", ], na2, assay(se3)["a", ]))
    expect_identical(assay(se)["b", ], c(assay(se1)["b", ], na2, assay(se3)["b", ]))
    expect_identical(assay(se)["c", ], c(assay(se1)["c", ], na2, na3))
    expect_identical(assay(se)["d", ], c(assay(se1)["d", ], na2, na3))
    expect_identical(assay(se)["e", ], c(assay(se1)["e", ], na2, na3))
    expect_identical(assay(se)["f", ], c(assay(se1)["f", ], na2, na3))
    expect_identical(assay(se)["h", ], c(assay(se1)["h", ], assay(se2)["h", ], na3))
    expect_identical(assay(se)["i", ], c(assay(se1)["i", ], assay(se2)["i", ], na3))
    expect_identical(assay(se)["j", ], c(assay(se1)["j", ], assay(se2)["j", ], assay(se3)["j", ]))
    expect_identical(assay(se)["k", ], c(na1, assay(se2)["k", ], assay(se3)["k", ]))
    expect_identical(assay(se)["l", ], c(na1, na2, assay(se3)["l", ]))
    expect_identical(assay(se)["m", ], c(na1, na2, assay(se3)["m", ]))
    expect_identical(assay(se)["n", ], c(na1, na2, assay(se3)["n", ]))
    ## test rowData
    expect_identical(rowData(se)$Prot, paste0("P", rownames(se)))
    expect_identical(sort(rowData(se)$Prot),
                     sort(union(c(rowData(se1)$Prot, rowData(se2)$Prot),
                                rowData(se3)$Prot)))
})

test_that("joinAssay", {
    jft <- joinAssays(feat2, 1:2)
    expect_identical(jft[["joinedAssay"]], QFeatures:::mergeSElist(list(se1, se2)))
    jft <- joinAssays(feat2, 2:1)
    expect_identical(jft[["joinedAssay"]], QFeatures:::mergeSElist(list(se2, se1)))
    jft <- joinAssays(feat2, c("assay1", "assay3"))
    expect_identical(jft[["joinedAssay"]], QFeatures:::mergeSElist(list(se1, se3)))
    jft <- joinAssays(feat2, 1:3)
    expect_identical(jft[["joinedAssay"]], QFeatures:::mergeSElist(list(se1, se2, se3)))
})

test_that("joinAssay errors", {
    expect_error(joinAssays(feat2, 1), "Need at least 2 assays to join")
    expect_error(joinAssays(feat2, "assay2"), "Need at least 2 assays to join")
    expect_error(joinAssay(se1, se2))
    expect_error(joinAssays(feat2, 1:2, name = "assay1"), "Assay with name 'assay1' already exists.")
    feat2 <- joinAssays(feat2, 1:3)
    expect_error(joinAssays(feat2, 1:3), "Assay with name 'joinedAssay' already exists.")
})


test_that("joinAssay return class", {
    m <- matrix(1:100, 10,
                dimnames = list(letters[1:10],
                                LETTERS[1:10]))
    library("SingleCellExperiment")
    sce <- SingleCellExperiment(list(m))
    se <- SummarizedExperiment(list(m))
    ## Joining two SCEs produces an SCE
    f <- QFeatures(ExperimentList(sce1 = sce[, 1:5], sce2 = sce[, 6:10]))
    expect_identical(class(joinAssays(f, i = 1:2)[["joinedAssay"]])[1],
                     "SingleCellExperiment")
    ## Joining two SEs produces an SE
    f <- QFeatures(ExperimentList(se1 = se[, 1:5], se2 = se[, 6:10]))
    expect_identical(class(joinAssays(f, i = 1:2)[["joinedAssay"]])[1],
                     "SummarizedExperiment")
    ## Joining an SE and an SCE throws an error
    f <- QFeatures(ExperimentList(se1 = se[, 1:5], sce2 = sce[, 6:10]))
    expect_error(joinAssays(f, i = 1:2))

})


test_that("aggregate and join (issue 81)", {
    ## See issue 81 for background on this unit test
    ## https://github.com/rformassspectrometry/QFeatures/issues/81
    data(feat2)
    feat2 <- aggregateFeatures(feat2, i = 1, name = "aggr1", fcol = "Prot", colSums)
    feat2 <- aggregateFeatures(feat2, i = 2, name = "aggr2", fcol = "Prot", colSums)
    feat2 <- aggregateFeatures(feat2, i = 3, name = "aggr3", fcol = "Prot", colSums)
    feat2 <- joinAssays(feat2, 1:3, "joinedAssay1")
    feat2 <- joinAssays(feat2, 4:6, "joinedAssay2")

    x <- feat2[["joinedAssay1"]]
    y <- feat2[["joinedAssay2"]]
    ## Both have same Prot variable
    expect_equal(rowData(x)[["Prot"]], rowData(y)[["Prot"]])
    ## rownames are different
    expect_equivalent(assay(x), assay(y))
})

test_that("join with an fcol", {
    ## See for background:
    ## https://github.com/rformassspectrometry/QFeatures/issues/221
    ## https://github.com/rformassspectrometry/QFeatures/pull/223
    
    ## Join assays that are not linked, 
    library(QFeatures)
    library(testthat)
    data(feat2)
    ## Create feature ids that overlap between sets
    rowData(feat2) <- lapply(rowData(feat2), function(x) {
        x$id <- seq_len(nrow(x))
        x
    })
    expect_message(
        test <- joinAssays(feat2, i = 1:3, fcol = "id"),
        regexp = "Using 'id' to join assays."
    )
    ## Check the joined assay has the expected dimensions, where
    ## number of rows should match the number of rows of largest set
    ## and the number of column should match the total number of
    ## colums
    coln <- unname(unlist(colnames(feat2)))
    rown <- as.character(seq_len(max(nrows(feat2))))
    ## Note that base::merge (used internally by joinAssays) does some
    ## weird reordering of the rows (although sort = FALSE)...
    ord <- c(1:7, 9:10, 8)
    expect_identical(
        dimnames(test[[4]]),
        list(rown[ord], coln)
    )
    ## Check the quant matrix is as expected
    exp_m <- matrix(
        NA, length(rown), length(coln), dimnames = list(rown, coln)
    )
    for (i in as.list(experiments(feat2))) {
        exp_m[rowData(i)$id, colnames(i)] <- assay(i)
    }
    exp_m <- exp_m[ord, ]
    expect_identical(exp_m, assay(test[[4]]))
    ## Check that the SE object is expected
    exp_cd <- do.call(rbind, lapply(experiments(feat2), function(x) {
        colData(x)[, "Var1", drop = FALSE]
    }))
    exp_se <- SummarizedExperiment(
        assays = list(exp_m), 
        rowData = DataFrame(id = as.integer(ord), row.names = ord),
        colData = exp_cd
    )
    expect_identical(exp_se, test[[4]])
    ## Check that the QF is expected
    exp_qf <- addAssay(feat2, exp_se, "joinedAssay")
    exp_qf <- addAssayLink(
        exp_qf, from = 1:3, 4, varFrom = rep("id", 3), varTo = "id"
    )
    expect_identical(exp_qf, test)
    ## When fcol has duplicated entries, the function automatically
    ## make the rownames unique and throws a warning
    rowData(feat2) <- lapply(rowData(feat2), function(x) {
        x$id <- c(1, seq_len(nrow(x) - 1))
        x
    })
    expect_message(
        expect_warning(
            test <- joinAssays(feat2, i = 1:3, fcol = "id"),
            regexp = "Duplicated entries found in ‘id’ in rowData of assay assay1; they are made unique."
        ),
        regexp = "Using 'id' to join assays."
    )
    ## Make sure this works when joining a subset of the assays
    data("feat3")
    expect_message(
        test <- joinAssays(feat3, i = 1:2, fcol = "Var"),
        regexp = "Using 'Var' to join assays."
    )
    exp_m <- assay(feat3[["psmsall"]])
    rownames(exp_m) <- sub("PSM", "", rownames(exp_m))
    expect_identical(assay(test[["joinedAssay"]]), exp_m)
    rd <- rowData(feat3[["psmsall"]])
    rownames(rd) <- sub("PSM", "", rownames(rd))
    exp_se <- SummarizedExperiment(
        assays = list(exp_m),
        rowData = rd,
        colData = colData(feat3[["psmsall"]])
    )
    expect_identical(test[["joinedAssay"]], exp_se)
    exp_qf <- addAssay(feat3, exp_se, "joinedAssay")
    exp_qf <- addAssayLink(
        exp_qf, from = 1:2, 8, varFrom = rep("Var", 2), varTo = "Var"
    )
    expect_identical(exp_qf, test)
})
