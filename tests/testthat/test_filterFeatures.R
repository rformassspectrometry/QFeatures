data(feat1)

test_that("filterFeatures", {
    feat1 <- aggregateFeatures(feat1, 1, fcol = "Sequence", name = "peptides", fun = colMedians)
    feat1 <- aggregateFeatures(feat1, 2, fcol = "Protein", name = "proteins", fun = colMedians)
    filter1 <- filterFeatures(feat1, ~  location == "Mitochondrion")
    filter2 <- filterFeatures(feat1, ~startsWith(location, "Mito"))
    filter3 <- filterFeatures(feat1, VariableFilter("location", "Mitochondrion"))
})
