data(feat1)
feat1 <- aggregateFeatures(feat1, "psms", "Sequence", "peptides")
feat2 <- addAssay(feat1[, , 1], feat1[["peptides"]], name = "peptides")

test_that(".createAssayLink", {
  expect_equivalent(feat1, 
                   .createAssayLink(feat2, from = "psms", to = "peptides",
                                    varFrom = "Sequence", varTo = "Sequence"))
  expect_equivalent(feat1, 
                    .createAssayLink(feat2, from = 1, to = 2,
                                     varFrom = "Sequence", varTo = "Sequence"))
  ## expect equivalent and not identical  because the `hits` in `feat1` are 
  ## sorted by `Sequence` by `aggregateFeatures`
})
