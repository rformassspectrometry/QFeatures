data("ft_na")

test_that("function: countUniqueFeatures", {
    ## Correct use
    ## No grouping
    expect_identical(countUniqueFeatures(ft_na,
                                         i = "na",
                                         colDataName = "counts")$counts,
                     c(2L, 3L, 4L))
    ## Grouping
    expect_identical(countUniqueFeatures(ft_na,
                                         i = "na",
                                         groupBy = "Y",
                                         colDataName = "counts")$counts,
                     c(2L, 2L, 2L))
    ## Expect errors
    ## colDataName already present in colData
    ft_na$foo <- 1
    expect_error(countUniqueFeatures(ft_na,
                                     i = "na",
                                     colDataName = "foo"),
                 regexp = "'foo' is already present in the colData.")
    ## The same sample is supplied twice
    ft_na <- addAssay(ft_na, ft_na[[1]][, 1, drop = FALSE], name = "na2")
    expect_error(countUniqueFeatures(ft_na,
                                     i = c("na", "na2")),
                 regexp = "The same sample is present in multiple assays.")
    
})