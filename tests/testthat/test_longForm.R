data(feat2)
colData(feat2)$X <- 1:12

test_that("longFormat is defunct", {
    expect_error(longFormat(feat2, rowvars = "Prot", colvars = "X"))
})

test_that("longForm", {
    ## Select a single colvars and rowvars
    lt <- longForm(feat2, rowvars = "Prot", colvars = "X")
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
    ## Select a single colvars and no rowvars (make sure that
    ## the implementation does not break the MAE implementation)
    lt <- longForm(feat2, colvars = "X")
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+1L)
    ## Select multiple rowvars
    lt <- longForm(feat2, rowvars = c("Prot", "x"))
    expect_equal(nrow(lt),
                 sum(apply(dims(feat2), 2, prod)))
    expect_identical(ncol(lt),
                     5L+2L)
    ## Test error: rowvars is missing in rowData
    expect_error(longForm(feat2, rowvars = "y"),
                 regexp = "not found")
})



longDF <- data.frame(
    rowname = rep(paste0("G", 1:10), 10),
    colname = rep(paste0("S", 1:10), each = 10),
    value = 1:100,
    assayName = "counts",
    colX = rep(LETTERS[1:10], each = 10),
    rowX = rep(letters[1:10], 10))
counts <- matrix(1:100, nrow = 10)
dimnames(counts) <- list(
    paste0("G", 1:10),
    paste0("S", 1:10))
cd <- data.frame(colX = LETTERS[1:10])
rd <- data.frame(rowX = letters[1:10])
se <- SummarizedExperiment(list(counts = counts),
                           colData = cd,
                           rowData = rd)
longDF <- as(longDF, "DataFrame")

test_that("longFormSE works with 1 assay", {
    expect_identical(longForm(se), longDF[, 1:4])
    expect_identical(longForm(se, colvars = "colX"), longDF[, 1:5])
    expect_identical(longForm(se, colvars = "colX", rowvars = "rowX"), longDF)
})


test_that("longFormSE works with 2 assays", {
    assay(se, "count2") <- assay(se)
    longDF2 <- longDF
    longDF2$assayName <- "count2"
    expect_identical(
        longForm(se, colvars = "colX", rowvars = "rowX"),
        rbind(longDF, longDF2))
})

test_that("longFormSE works with unnamed assays", {
    ## one "" names
    assay(se, 2) <- assay(se)
    longDF2 <- longDF
    longDF2$assayName <- ""
    expect_identical(
        longForm(se, colvars = "colX", rowvars = "rowX"),
        rbind(longDF, longDF2))
    alst <- assays(se)
    names(alst) <- NULL
    assays(se) <- alst
    longDF$assayName <- 1L
    longDF2$assayName <- 2L
    ## No names
    expect_equal(
        longForm(se, colvars = "colX", rowvars = "rowX"),
        rbind(longDF, longDF2))
})