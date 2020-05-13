## Define numeric matrices
m1 <- matrix(1:40, ncol = 4) + 0.0
m2 <- matrix(1:16, ncol = 4) + 0.0
m3 <- matrix(1:28, ncol = 4) + 0.0
colnames(m1) <- paste0("S", 1:4)
colnames(m2) <- paste0("S", 5:8)
colnames(m3) <- paste0("S", 9:12)
rownames(m1) <- letters[1:10]
rownames(m2) <- letters[8:11]
rownames(m3) <- letters[c(1, 2, 10:14)]

df1 <- DataFrame(Prot = paste0("P", rownames(m1)),
                 x = rnorm(1:10),
                 row.names = rownames(m1))
cd1 <- DataFrame(Var1 = rnorm(4),
                 Var2 = LETTERS[1:4],
                 row.names = colnames(m1))


df2 <- DataFrame(Prot = paste0("P", rownames(m2)),
                 x = rnorm(1:4),
                 y = rnorm(1:4),
                 row.names = rownames(m2))
cd2 <- DataFrame(Var1 = rnorm(4), 
                 Var2 = LETTERS[5:8],
                 row.names = colnames(m2))

df3 <- DataFrame(Prot = paste0("P", rownames(m3)),
                 x = rnorm(1:7),
                 y = rnorm(1:7),
                 row.names = rownames(m3))
cd3 <- DataFrame(Var1 = rnorm(4),
                 row.names = colnames(m3))


se1 <- SummarizedExperiment(m1, df1, colData = cd1)
se2 <- SummarizedExperiment(m2, df2, colData = cd2)
se3 <- SummarizedExperiment(m3, df3, colData = cd3)

el <- list(assay1 = se1, assay2 = se2, assay3 = se3)
ft <- Features(el)


test_that("merge SEs (1,2)", {
    se <- mergeSElist(list(se1, se2))
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
    se <- mergeSElist(list(se2, se3))
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
    se <- mergeSElist(list(se3, se2))
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
    se <- mergeSElist(list(se1, se2, se3))
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
    jft <- joinAssays(ft, 1:2)
    expect_identical(jft[["joinedAssay"]], mergeSElist(list(se1, se2)))
    jft <- joinAssays(ft, 2:1)
    expect_identical(jft[["joinedAssay"]], mergeSElist(list(se2, se1)))
    jft <- joinAssays(ft, c("assay1", "assay3"))
    expect_identical(jft[["joinedAssay"]], mergeSElist(list(se1, se3)))
    jft <- joinAssays(ft, 1:3)
    expect_identical(jft[["joinedAssay"]], mergeSElist(list(se1, se2, se3)))    
})
