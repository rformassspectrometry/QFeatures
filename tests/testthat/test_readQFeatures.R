data(hlpsms)
x <- hlpsms
x$file <- paste0("File", sample(1:3, nrow(x), replace = TRUE))
colann <- data.frame(file = rep(paste0("File", 1:3), each = 10),
                     Channel = rep(names(x)[1:10], 3))


test_that(".checkWarnEcol works", {
    expect_identical(.checkWarnEcol(NULL, NULL), NULL)
    expect_identical(.checkWarnEcol(1, NULL), 1)
    expect_warning(res <- .checkWarnEcol(NULL, 1),
                   "'ecol' is deprecated, use 'quantCols' instead")
    expect_error(.checkWarnEcol(1, 1),
                 "Use 'quantCols' only.")
    expect_identical(res, 1)
})

test_that(".checkQuantCols works", {
    expect_error(readQFeatures(x),
                 "Provide one of 'colData' or 'quantCols'")
    expect_error(readQFeatures(x, colData = data.frame()),
                 "'colData' must contain a column called 'quantCols'")
    ## expect_warning(readQFeatures(x, colData = data.frame(quantCols = 1:10),
    ##                              quantCols = 1:10),
    ##                "Ignoring 'quantCols', using 'colData'.")
    expect_error(readQFeatures(x, colData = data.frame(quantCols = "foo")),
                 "Some column names in 'quantCols' are not found in 'assayData'")
    expect_error(readQFeatures(x, colData = data.frame(quantCols = c("X126", "foo"))),
                 "Some column names in 'quantCols' are not found in 'assayData'")
})

test_that(".checkRunCol works", {
    ## Single-set case
    expect_null(.checkRunCol(, , NULL))
    ##################################################
    ## Multi-set case and runCol: positive control
    expect_identical(.checkRunCol(x, NULL, "file"), x[["file"]])
    ## --------------------------------------------
    ## Multi-set case and runCol: errors
    expect_error(.checkRunCol(x, NULL, c("file", "file")),
                 "'runCol' must contain the name of a single column in 'assayData'.")
    expect_error(.checkRunCol(x, NULL, "nofile"),
                 "not found in 'assayData'.")
    ##################################################
    ## Multi-set case, runCol and colData: positive control
    colann <- data.frame(runCol = rep(paste0("File", 1:3), each = 10),
                         quantCol = rep(names(x)[1:10], 3))
    expect_identical(.checkRunCol(x, colann, "file"), x[["file"]])
    ## --------------------------------------------
    ## Multi-set case, runCol and colData: errors
    expect_error(.checkRunCol(x, colann, "nofile"),
                 "not found in 'assayData'.")
    names(colann)[1] <- "noRunCol"
    expect_error(.checkRunCol(x, colann, "file"),
                 "'colData' must contain a column called 'runCol'.")
    colann <- data.frame(runCol = rep(paste0("File", 1:2), each = 15),
                         quantCol = rep(names(x)[1:10], 3))
    expect_warning(.checkRunCol(x, colann, "file"),
                   "Some runs are missing in 'colData':")
})

test_that("readQFeatures: colData and quantCols are equivalent", {
    r1 <- readQFeatures(x, colData = data.frame(quantCols = names(x)[1:10]))
    r2 <- readQFeatures(x, quantCols = 1:10)
    colData(r1) <- colData(r2) ## ignore colData
    expect_identical(r1, r2)
})


test_that("readQFeatures: testing use cases", {
    require("SummarizedExperiment")
    require("MultiAssayExperiment")
    ## cf this comment for an overview of the use 3 cases:
    ## https://github.com/rformassspectrometry/QFeatures/pull/200#issuecomment-1986862723

    ## Case 1: Single-set, multiple quantitative cols
    se_exp <- readSummarizedExperiment(x, 1:10)
    colAnnot <- DataFrame(
        quantCols = names(x)[1:10],
        file = rep(paste0("File", 1:3), length.out = 10),
        row.names = grep("^X", colnames(x), value = TRUE)
    )
    ## Without colAnnot
    expect_identical(
        readQFeatures(x, quantCols = 1:10),
        QFeatures(List(quants = se_exp),
                  metadata = list("._type" = "bulk"))
    )
    ## With colAnnot
    expect_identical(
        readQFeatures(x, 1:10, colData = colAnnot),
        QFeatures(List(quants = se_exp), colData = colAnnot,
                  metadata = list("._type" = "bulk"))
    )
    ## Without quantCols
    expect_identical(
        qf <- readQFeatures(x, colData = colAnnot),
        QFeatures(List(quants = se_exp), colData = colAnnot,
                  metadata = list("._type" = "bulk"))
    )
    ## Check colnames
    expect_identical(
        colnames(qf),
        CharacterList(quants = colnames(se_exp))
    )
    ## With colAnnot and quantCols, but with different x columns and quantCols (PR #234)
    shuffledSamplesOrder <- sample(1:10, 10L)
    shuffledColAnnot <- colAnnot[shuffledSamplesOrder, , drop = FALSE]
    expect_identical(
        readQFeatures(x, shuffledSamplesOrder, colData = shuffledColAnnot),
        QFeatures(List(quants = se_exp), colData = colAnnot,
                  metadata = list("._type" = "bulk"))
    )
    ## With colAnnot and quantCols, but with different colAnnot row order (PR #234)
    expect_identical(
        readQFeatures(x, 1:10, colData = shuffledColAnnot),
        QFeatures(List(quants = se_exp), colData = colAnnot,
                  metadata = list("._type" = "bulk"))
    )

    ## Case 2: Multiple-set, one quantitative col
    se_exp <- readSummarizedExperiment(x, 1)
    colAnnot <- DataFrame(
        file = unique(x$file),
        row.names = unique(x$file)
    )
    el_exp <- lapply(paste0("File", 1:3), function(file) {
        out <- se_exp[rowData(se_exp)$file == file, ]
        colnames(out) <- file
        out
    })
    names(el_exp) <- paste0("File", 1:3)
    ## Without colAnnot
    expect_identical(
        readQFeatures(x, quantCol = 1, runCol = "file"),
        QFeatures(List(el_exp),
                  metadata = list("._type" = "bulk"))
    )
    ## With colAnnot
    colAnnot$runCol <- colAnnot$file
    expect_identical(
        readQFeatures(x, quantCols = 1, runCol = "file", colData = colAnnot),
        QFeatures(List(el_exp), colData = colAnnot[paste0("File", 1:3), ],
                  metadata = list("._type" = "bulk"))
    )
    ## When length(quantCol) == 1, colAnnot doesn't need a quantCol
    expect_identical(
        readQFeatures(x, 1, runCol = "file", colData = colAnnot),
        QFeatures(List(el_exp), colData = colAnnot[paste0("File", 1:3), ],
                  metadata = list("._type" = "bulk"))
    )
    ## Without quantCols
    colAnnot$quantCols <- colnames(x)[1]
    expect_identical(
        qf <- readQFeatures(x, runCol = "file", colData = colAnnot),
        QFeatures(List(el_exp), colData = colAnnot[paste0("File", 1:3), ],
                  metadata = list("._type" = "bulk"))
    )
    ## Check colnames
    expect_identical(
        colnames(qf),
        CharacterList(File1 = "File1", File2 = "File2", File3 = "File3")
    )

    ## Case 3: Multiple-set, multiple quantitative cols
    se_exp <- readSummarizedExperiment(x, 1:10)
    runs <- rep(unique(x$file), each = ncol(se_exp))
    quantCols <- rep(colnames(se_exp), length(unique(x$file)))
    colAnnot <- DataFrame(
        file = runs,
        quantCols = quantCols,
        row.names = paste0(runs, "_", quantCols)
    )
    el_exp <- lapply(paste0("File", 1:3), function(file) {
        out <- se_exp[rowData(se_exp)$file == file, ]
        colnames(out) <- paste0(file, "_", colnames(out))
        out
    })
    names(el_exp) <- paste0("File", 1:3)
    ## Without colAnnot
    expect_identical(
        readQFeatures(x, quantCols = 1:10, runCol = "file"),
        QFeatures(List(el_exp),
                  metadata = list("._type" = "bulk"))
    )
    ## With colAnnot
    colAnnot$runCol <- colAnnot$file
    expect_identical(
        readQFeatures(x, quantCols = 1:10, runCol = "file", colData = colAnnot),
        QFeatures(List(el_exp), colData = colAnnot[order(rownames(colAnnot)), ],
                  metadata = list("._type" = "bulk"))
    )
    ## Without quantCols
    expect_identical(
        qf <- readQFeatures(x, runCol = "file", colData = colAnnot),
        QFeatures(List(el_exp), colData = colAnnot[order(rownames(colAnnot)), ],
                  metadata = list("._type" = "bulk"))
    )
    ## Check colnames
    expect_identical(
        colnames(qf),
        CharacterList(split(rownames(colAnnot), colAnnot$file))
    )

    ## Test removeEmptyCols
    x_na <- x
    x_na$X126 <- NA
    se_na <- readSummarizedExperiment(x_na, 1:10)
    expect_identical(
        readQFeatures(x_na, quantCol = 1:10, removeEmptyCols = FALSE),
        QFeatures(List(quants = se_na),
                  metadata = list("._type" = "bulk"))
    )
    expect_identical(
        readQFeatures(x_na, quantCols = 1:10, removeEmptyCols = TRUE),
        QFeatures(List(quants = se_na[, -1]),
                  metadata = list("._type" = "bulk"))
    )
})

test_that("readQFeatures: test polymorphism", {
    se_exp <- readSummarizedExperiment(x, 1:10)
    ## Test name argument
    expect_identical(
        readQFeatures(x, quantCols = 1:10), ## no name = default to quants
        readQFeatures(x, quantCols = 1:10, name = "quants")
    )
    expect_identical(
        readQFeatures(x, quantCols = 1:10, name = "foo"),
        QFeatures(List(foo = se_exp),
                  metadata = list("._type" = "bulk"))
    )
    ## Test quantCols polymorphism
    expect_identical(
        readQFeatures(x, quantCols = colnames(x)[1:10]),
        readQFeatures(x, quantCols = seq_len(ncol(x)) %in% 1:10)
    )
    ## colAnnot is coercible to data.frame
    colAnnot <- DataFrame(
        quantCols = colnames(se_exp),
        file = rep(paste0("File", 1:3), length.out = ncol(se_exp)),
        row.names = colnames(se_exp)
    )
    expect_identical(
        readQFeatures(x, colData = colAnnot),
        readQFeatures(x, colData = as.list(colAnnot))
    )
})


test_that("readQFeatures: errors, warnings and messages", {
    se_exp <- readSummarizedExperiment(x, 1:10)
    colAnnot <- DataFrame(
        quantCols = names(x)[1:10],
        file = rep(paste0("File", 1:3), length.out = 10),
        row.names = grep("^X", colnames(x), value = TRUE)
    )
    ## no quantCols and no colAnnot = error
    expect_error(
        readQFeatures(x),
        regexp = "Provide one of 'colData' or 'quantCols', both mustn't be NULL."
    )
    ## some quantCols are missing = error
    expect_error(
        readQFeatures(x, quantCols = "foo"),
        regexp = "Some column names in 'quantCols' are not found in 'assayData': foo"
    )
    ## if colAnnot is not NULL, it must contain a quantCols column
    colAnnot$quantCols <- NULL
    expect_error(
        readQFeatures(x, quantCols = 1:10, colData = colAnnot),
        regexp = "^'colData' must contain a column called 'quantCols'"
    )
    expect_error(
        readQFeatures(x, colData = colAnnot),
        regexp = "When 'quantCols' is NULL, 'colData' must contain a column called 'quantCols'"
    )
    expect_error(
        readQFeatures(x, quantCols = 1:10, runCol = "file", colData = colAnnot),
        regexp = "^'colData' must contain a column called 'quantCols'"
    )
    ## runCol is a vector = error
    expect_error(
        readQFeatures(x, quantCols = 1:10, runCol = x$file),
        regexp = "'runCol' must contain the name of a single column in 'assayData'."
    )
    ## runCol not found in assay data = error
    expect_error(
        readQFeatures(x, quantCols = 1:10, runCol = "foo"),
        regexp = "'foo' .provided as 'runCol'. not found in 'assayData'."
    )
    ## runCol missing in colAnnot = error
    colAnnot <- DataFrame(
        quantCols = names(x)[1:10],
        file = rep(paste0("File", 1:3), length.out = 10),
        row.names = grep("^X", colnames(x), value = TRUE)
    )
    expect_error(
        readQFeatures(x, quantCols = 1:10, runCol = "file", colData = colAnnot),
        regexp = "When 'runCol' is not NULL, 'colData' must contain a column called 'runCol'."
    )

    ## Missing runs in colAnnot = warning
    se_exp <- readSummarizedExperiment(x, 1)
    colAnnot <- DataFrame(
        runCol = paste0("File", 1:2),
        row.names = paste0("File", 1:2)
    )
    el_exp <- lapply(paste0("File", 1:3), function(file) {
        out <- se_exp[rowData(se_exp)$file == file, ]
        colnames(out) <- file
        out
    })
    names(el_exp) <- paste0("File", 1:3)
    expect_warning(
        qf <- readQFeatures(x, quantCols = 1, runCol = "file", colData = colAnnot),
        regexp = "Some runs are missing in 'colData': File3"
    )
    ## Missing annotations are autatomically filled with NA
    colAnnot["File3", ] <- NA
    expect_identical(qf, QFeatures(List(el_exp), colData = colAnnot,
                                   metadata = list("._type" = "bulk")))

    ## Test verbose = messages
    expect_no_message(readQFeatures(x, quantCols = 1:10, verbose = FALSE))
    msgs <- capture_messages(readQFeatures(x, quantCols = 1:10))
    expect_identical(msgs, c(
        "Checking arguments.\n",
        "Loading data as a 'SummarizedExperiment' object.\n",
        "Formatting sample annotations (colData).\n",
        "Formatting data as a 'QFeatures' object.\n"
    ))
    msgs <- capture_messages(readQFeatures(x, quantCols = 1:10, runCol = "file"))
    expect_identical(msgs, c(
        "Checking arguments.\n",
        "Loading data as a 'SummarizedExperiment' object.\n",
        "Splitting data in runs.\n",
        "Formatting sample annotations (colData).\n",
        "Formatting data as a 'QFeatures' object.\n"
    ))
})

test_that("readSummarizedExperiment", {
    expect_warning(
        ft1 <- readSummarizedExperiment(x, ecol = 1:10),
        regexp = "'ecol' is deprecated, use 'quantCols' instead."
    )
    ft2 <- readSummarizedExperiment(x, quantCols = 1:10)
    expect_equal(ft1, ft2)
    expect_warning(
        ft3 <- readSummarizedExperiment(x, ecol = 1:10, fnames = "Sequence"),
        regexp = "'ecol' is deprecated, use 'quantCols' instead."
    )
    ft4 <- readSummarizedExperiment(x, quantCols = 1:10, fnames = "Sequence")
    expect_equal(ft3, ft4)
    ft5 <- readSummarizedExperiment(x, quantCols = 1:10, fnames = 11)
    expect_equal(ft3, ft5)
    ## Read data with only 1 quantitation column
    ft5 <- readSummarizedExperiment(x, quantCols = 1, fnames = 11)
    ## Check column names
    quantCols <- c("X126", "X127C", "X127N", "X128C", "X128N", "X129C",
                   "X129N", "X130C", "X130N", "X131")
    expect_identical(colnames(ft3), quantCols)
    expect_identical(colnames(ft5), quantCols[1])
    ## Provide quantCols as logical
    quantCols <- seq_along(x) %in% 1:10
    expect_identical(ft1, readSummarizedExperiment(x, quantCols = quantCols))
    ## Expect errors
    quantCols <- LETTERS[1:10]
    expect_error(
        readSummarizedExperiment(x, quantCols = quantCols, name = "psms"),
        regexp = "Column identifiers A, B, C, D, E, F, G, H, I, J not recognised among"
    )
    expect_error(
        readSummarizedExperiment(x, quantCols = 1:10, fnames = "not_present"),
        regexp = "not_present not found among"
        )
    expect_true(inherits(ft1, "SummarizedExperiment"))
})

test_that(".splitSE", {
    m <- matrix(1:100, ncol = 10,
                dimnames = list(paste0("row", 1:10),
                                paste0("col", 1:10)))
    se <- SummarizedExperiment(assay = m, 
                               rowData = DataFrame(rowDataCol = 1:nrow(m)%%3),
                               colData = DataFrame(colvar = 1:ncol(m)%%5))
    ## Split by row
    expect_identical(length(.splitSE(se, "rowDataCol")), 3L)
    ## Split by col
    expect_identical(length(.splitSE(se, "colvar")), 5L)
    ## Error: variable not found
    expect_error(.splitSE(se, "foo"),
                 regexp = "not found")
    ## Error: factor is too short
    expect_error(.splitSE(se, factor(1:3)),
                 regexp = "not compatible with dim")
})
