data(feat2)
colData(feat2)$X <- 1:12

test_that("longFormat is deprecated", {
    expect_warning(longFormat(feat2, rowvars = "Prot", colvars = "X"),
                   "'longFormat' is deprecated.")
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

## longDF <- data.frame(
##     rowname = rep(paste0("G", 1:10), 10),
##     colname = rep(paste0("S", 1:10), each = 10),
##     value = 1:100,
##     assayName = "counts",
##     colX = rep(LETTERS[1:10], each = 10),
##     rowX = rep(letters[1:10], 10))

## counts <- matrix(1:100, nrow = 10)
## dimnames(counts) <- list(
##     paste0("G", 1:10),
##     paste0("S", 1:10))
## cd <- data.frame(colX = LETTERS[1:10])
## rd <- data.frame(rowX = letters[1:10])
## se <- SummarizedExperiment(list(counts = counts),
##                            colData = cd,
##                            rowData = rd)


## all.equal(longFormSE(se), longDF[, 1:4])
## all.equal(longFormSE(se, colvars = "colX"), longDF[, 1:5])
## all.equal(longFormSE(se, colvars = "colX", rowvars = "rowX"), longDF)