data(feat1)

test_that("filterFeatures", {
    feat1 <- aggregateFeatures(feat1, 1, fcol = "Sequence", name = "peptides", fun = colMedians)
    feat1 <- aggregateFeatures(feat1, 2, fcol = "Protein", name = "proteins", fun = colMedians)
    ## Test character filters
    filter1 <- filterFeatures(feat1, ~  location == "Mitochondrion")
    filter2 <- filterFeatures(feat1, ~startsWith(location, "Mito"))
    filter3 <- filterFeatures(feat1, VariableFilter("location", "Mitochondrion"))
    filter4 <- filterFeatures(feat1, VariableFilter("location", "unknown", condition = "!="))
    filter5 <- filterFeatures(feat1, VariableFilter("location", "unknown", condition = "==", not = TRUE))
    filter6 <- filterFeatures(feat1, ~ location != "unknown")
    filter7 <- filterFeatures(feat1, VariableFilter("location", "ochon", condition = "contains"))
    expect_equal(filter1, filter2)
    expect_equal(filter1, filter3)
    expect_equal(filter1, filter4)
    expect_equal(filter1, filter5)
    expect_equal(filter1, filter6)
    expect_equal(filter1, filter7)
    expect_identical(lengths(filter1), c(6L, 2L, 1L))
    ## Test numerical filters
    filter1 <- filterFeatures(feat1, VariableFilter("pval", 0.03, "<="))
    filter2 <- filterFeatures(feat1, ~ pval <= 0.03)
    expect_equal(filter1, filter2)
    ## Test no match filters
    filter1 <- filterFeatures(feat1, ~  location != "Mitochondrion")
    expect_identical(lengths(filter1), c(4L, 1L, 1L))
    expect_true(isEmpty(filterFeatures(feat1, VariableFilter("location", "not"))))
    expect_true(isEmpty(filterFeatures(feat1, ~ location == "not")))
    expect_true(isEmpty(filterFeatures(feat1, VariableFilter("foo", "bar"))))
    expect_true(isEmpty(filterFeatures(feat1, ~ foo == "bar")))
    expect_true(isEmpty(filterFeatures(feat1, ~ is.na(pval))))
    ## Test fraud filters
    expect_error(VariableFilter("pval", TRUE, "<="))
    expect_error(VariableFilter("location", TRUE, "!="))
})  

