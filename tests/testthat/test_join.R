data(feat2)

se1 <- feat2[[1]]
se2 <- feat2[[2]]
se3 <- feat2[[3]]


test_that("merge SEs (1,2)", {
    se <- Features:::mergeSElist(list(se1, se2))
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
    se <- Features:::mergeSElist(list(se2, se3))
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
    se <- Features:::mergeSElist(list(se3, se2))
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
    se <- Features:::mergeSElist(list(se1, se2, se3))
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
    expect_identical(jft[["joinedAssay"]], Features:::mergeSElist(list(se1, se2)))
    jft <- joinAssays(feat2, 2:1)
    expect_identical(jft[["joinedAssay"]], Features:::mergeSElist(list(se2, se1)))
    jft <- joinAssays(feat2, c("assay1", "assay3"))
    expect_identical(jft[["joinedAssay"]], Features:::mergeSElist(list(se1, se3)))
    jft <- joinAssays(feat2, 1:3)
    expect_identical(jft[["joinedAssay"]], Features:::mergeSElist(list(se1, se2, se3)))    
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
    f <- Features(ExperimentList(sce1 = sce[, 1:5], sce2 = sce[, 6:10]))
    expect_identical(class(joinAssays(f, i = 1:2)[["joinedAssay"]])[1],
                     "SingleCellExperiment")
    ## Joining two SEs produces an SE
    f <- Features(ExperimentList(se1 = se[, 1:5], se2 = se[, 6:10]))
    expect_identical(class(joinAssays(f, i = 1:2)[["joinedAssay"]])[1],
                     "SummarizedExperiment")
    ## Joining an SE and an SCE throws an error
    f <- Features(ExperimentList(se1 = se[, 1:5], sce2 = sce[, 6:10]))
    expect_error(joinAssays(f, i = 1:2))

})
