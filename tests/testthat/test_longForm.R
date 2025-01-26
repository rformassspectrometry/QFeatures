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
