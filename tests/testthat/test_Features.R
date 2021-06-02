data(feat1)

test_that("empty QFeatures", {
    feat0 <- QFeatures()
    expect_true(validObject(feat0))
    expect_true(isEmpty(feat0))
    expect_null(show(feat0))
})

test_that("Manual QFeatures", {
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
    feat2 <- QFeatures(list(psms = psms), colData = coldata)
    expect_true(validObject(feat2))
    ## subsetting
    expect_null(show(feat2))
    expect_equal(psms, feat2[[1]])
    expect_equal(psms, feat2[["psms"]])
    expect_equal(feat2, feat2[1:10, 1:2, 1])
    ## compare to serialised data
    ## data(feat1)
    ## expect_equivalent(feat1, feat2)
})


test_that("addAssay", {
    data(feat1)
    assay2 <- feat1[[1]]
    feat1 <- addAssay(feat1, assay2, name = "psms2")
    expect_identical(names(feat1), c("psms", "psms2"))
    expect_identical(feat1[[1]], feat1[[2]])
})

test_that("[,QFeatures", {
    data(feat1)
    feat1 <- aggregateFeatures(feat1, "psms", "Sequence", "peptides")
    expect_true(expect_warning(validObject(feat1[, , "psms"]),
                               regexp = "'experiments' dropped; see 'metadata'"))
    expect_true(expect_warning(validObject(feat1[, , "peptides"]),
                               regexp = "'experiments' dropped; see 'metadata'"))
    expect_true(expect_warning(validObject(feat1[, , 1]),
                               regexp = "'experiments' dropped; see 'metadata'"))
    expect_true(expect_warning(validObject(feat1[, , 2]),
                               regexp = "'experiments' dropped; see 'metadata'"))
})


test_that("RowData", {
    data(feat2)
    rd <- rowData(feat2)
    expect_identical(names(rd), paste0("assay", 1:3))
    expect_identical(rd[[1]], rowData(feat2[[1]]))
    expect_identical(rd[[2]], rowData(feat2[[2]]))
    expect_identical(rd[[3]], rowData(feat2[[3]]))
})


test_that("rowDataNames", {
    rdn <- rowDataNames(feat1)
    expect_identical(length(feat1), length(rdn))
    expect_identical(names(feat1), names(rdn))
    for (i in seq_along(length(feat1)))
        expect_identical(rdn[[i]], names(rowData(feat1[[i]])))
})


test_that("selectRowData", {
    x <- c("Sequence", "Protein")
    ft <- selectRowData(feat1, x)
    expect_identical(length(ft), length(feat1))
    expect_identical(names(ft), names(feat1))
    expect_identical(rowDataNames(ft)[[1]], x)
    expect_error(selectRowData(feat1))
    expect_message(ft <- selectRowData(feat1,
                                       c("Sequence", "Protein",
                                         "Var", "location", "pval",
                                         "var_not_found")))
    expect_identical(feat1, ft)
})


test_that("getRowData", {
    data(feat2)
    ## Get 1 variable from 1 assay
    rd <- getRowData(feat2, 1, "Prot")
    expect_true(inherits(rd, "DFrame"))
    expect_identical(unique(rd$assay), names(feat2)[1])
    expect_identical(colnames(rd), c("assay", "rowname", "Prot"))
    expect_identical(nrow(rd), sum(dims(feat2)[1, 1]))
    ## Get several variables from 1 assay
    rd <- getRowData(feat2, 1, c("Prot", "x"))
    expect_true(inherits(rd, "DFrame"))
    expect_identical(unique(rd$assay), names(feat2)[1])
    expect_identical(colnames(rd), c("assay", "rowname", "Prot", "x"))
    expect_identical(nrow(rd), sum(dims(feat2)[1, 1]))
    ## Get 1 variable from several assays
    rd <- getRowData(feat2, 1:2, "Prot")
    expect_true(inherits(rd, "DFrame"))
    expect_identical(unique(rd$assay), names(feat2)[1:2])
    expect_identical(colnames(rd), c("assay", "rowname", "Prot"))
    expect_identical(nrow(rd), sum(dims(feat2)[1, 1:2]))
    ## Get several variables from several assays
    rd <- getRowData(feat2, 1:2, c("Prot", "x"))
    expect_true(inherits(rd, "DFrame"))
    expect_identical(unique(rd$assay), names(feat2)[1:2])
    expect_identical(colnames(rd), c("assay", "rowname", "Prot", "x"))
    expect_identical(nrow(rd), sum(dims(feat2)[1, 1:2]))
    ## Get all common variable from all assays (missing rowDataCols)
    rd <- getRowData(feat2, i = seq_along(feat2))
    expect_true(inherits(rd, "DFrame"))
    expect_identical(unique(rd$assay), names(feat2))
    expect_identical(colnames(rd), c("assay", "rowname", "Prot", "x"))
    expect_identical(nrow(rd), sum(dims(feat2)[1, ]))
    ## Error: variable not present in one of the assays 
    expect_error(getRowData(feat2, seq_along(feat2), "y"),
                 regexp = "Some 'rowDataCols' are not found in the rowData")
})


test_that("setRowData", {
    data(feat2)
    ## Add rowData
    ## In one assay
    repl1 <- lapply(rowData(feat2), function(x) DataFrame(foo = rep("bar", nrow(x))))
    feat3 <- setRowData(feat2, rowDataCols = "foo", replacement = repl1[1])
    expect_identical(unique(getRowData(feat3, 1, "foo")$foo), "bar")
    ## In more assay
    feat3 <- setRowData(feat2, rowDataCols = "foo", replacement = repl1[1:2])
    expect_identical(unique(getRowData(feat3, names(repl1[1:2]), "foo")$foo), "bar")
    expect_identical(unique(getRowData(feat3, names(repl1[1:2]), "foo")$assay), names(repl1[1:2]))
    ## More variables
    repl2 <- lapply(repl1, function(x) cbind(x, foo2 = x$foo))
    feat3 <- setRowData(feat2, rowDataCols = c("foo", "foo2"), replacement = repl2[1:2])
    expect_true(all(any(rowDataNames(feat3)[1:2] == "foo2")))
    ## Update rowData
    repl3 <- lapply(repl2, function(x){
        x$foo2 <- "bar2"
        x
    })
    feat3 <- setRowData(feat3, rowDataCols = "foo2", replacement = repl3[1:2])
    expect_identical(unique(getRowData(feat3, names(repl3[1:2]), "foo2")$foo2), "bar2")
    ## Remove rowData
    ## Remove multiple columns existing in some assay, but absent in others
    repl4 <- lapply(rowData(feat2), function(x) NULL)
    feat4 <- setRowData(feat3, rowDataCols = c("foo", "foo2"), replacement = repl4)
    expect_identical(rowData(feat2), rowData(feat4))
    ## Remove in one assay and add in another
    repl4[[3]] <- repl3[[3]]
    feat4 <- setRowData(feat3, rowDataCols = c("foo", "foo2"), replacement = repl4)
    expect_identical(rowData(feat2)[1:2], rowData(feat4)[1:2])
    expect_identical(unique(getRowData(feat4, names(repl4[3]), "foo2")$foo2), "bar2")
    ## Error when replacement names point to an assay that is not in 
    ## the QFeatures object
    repl5 <- repl4
    names(repl5) <- paste0("foo", length(repl5))
    expect_error(setRowData(feat2, rowDataCols = c("foo", "foo2"), replacement = repl5),
                 regexp = "replacement.*in.*object.* is not TRUE")
    ## Error rowDataCols not in column names of replacement tables
    expect_error(setRowData(feat2, rowDataCols = c("foo3"), replacement = repl3[1]),
                 regexp = "invalid names")
    ## Warning replacement table with wrong dimensions
    repl6 <- repl3
    repl6[[1]] <- repl6[[1]][1:7, ]
    expect_warning(setRowData(feat2, rowDataCols = c("foo"), replacement = repl6[1]),
                 regexp = "number of values to be replaced")
})


test_that("renaming", {
    data(feat1)
    feat1 <- aggregateFeatures(feat1, "psms", fcol = "Sequence",
                               name = "peptides", fun = colSums)
    feat1 <- aggregateFeatures(feat1, "peptides", fcol = "Protein",
                               name = "proteins", fun = colSums)
    expect_true(validObject(feat1))
    feat2 <- feat1
    names(feat2) <- LETTERS[1:3]
    expect_true(validObject(feat2))
    # Expect errors
    feat2@assayLinks[[1]]@name <- "foo"
    expect_error(validObject(feat2),
                 regexp = "@names not valid")
    feat2 <- feat1
    feat2@assayLinks[[1]]@from <- "bar"
    expect_error(validObject(feat2),
                 regexp = "@from not valid")
    expect_error(names(feat2) <- 1:3,
                 regexp = "must be a character")
    expect_error(names(feat2) <- letters[c(1,2,2)],
                 regexp = "is duplicated$")
})

test_that("assays must have unique rownames", {
    hlpsms <- hlpsms[1:10, ]
    ft1 <- readQFeatures(hlpsms, ecol = 1:10, name = "psms", fname = "Sequence")
    rownames(ft1[[1]][1:2]) <- rep("1", 2)
    expect_error(validObject(ft1))
})


test_that("longFormat", {
    colData(feat2)$X <- 1:12
    ## Select a single colDataCols and rowDataCols
    lt <- longFormat(feat2, rowDataCols = "Prot", colDataCols = "X")
    ## Check dimensions
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+2L)
    ## Check content
    expect_identical(lt$Prot,
                     unname(unlist(lapply(rowData(feat2), 
                                          ## Repeat 4x because 4 samples
                                          function(x) rep(x$Prot, 4)))))
    ## Select a single colDataCols and no rowDataCols (make sure that 
    ## the implementation does not break the MAE implementation)
    lt <- longFormat(feat2, colDataCols = "X")
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+1L)
    ## Select multiple rowDataCols
    lt <- longFormat(feat2, rowDataCols = c("Prot", "x"))
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+2L)
    ## Test error: rowDataCols is missing in rowData
    expect_error(longFormat(feat2, rowDataCols = "y"),
                 regexp = "not found")
})

