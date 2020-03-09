data(feat1)
psms <- feat1[["psms"]]
feat1 <- aggregateFeatures(feat1, "psms", "Sequence", "peptides")
feat2 <- addAssay(feat1[, , 1], feat1[["peptides"]], name = "peptides")


test_that(".createAssayLink", {
  ## Correct use
  ## we here use expect_equivalent and not identical  because the `hits` in 
  ## `feat1` are sorted by `Sequence` by `aggregateFeatures`
  expect_equivalent(feat1@assayLinks[[2]], 
                   .createAssayLink(feat2, from = "psms", to = "peptides",
                                    varFrom = "Sequence", varTo = "Sequence"))
  
  ## No match between feature variables
  expect_error(.createAssayLink(feat2, from = "psms", to = "peptides",
                                varFrom = "Sequence", varTo = "Protein"))
})

test_that("createAssayLink", {
  expect_equivalent(feat1, 
                    createAssayLink(feat2, from = "psms", to = "peptides",
                                     varFrom = "Sequence", varTo = "Sequence"))
  expect_equivalent(feat1, 
                    createAssayLink(feat2, from = 1, to = 2,
                                     varFrom = "Sequence", varTo = "Sequence"))
  expect_equivalent(feat1, 
                    createAssayLink(feat2, from = 1, to = 2, 
                                    varFrom = "Sequence"))
})


test_that("createOneToOneAssayLink", {
  feat3 <- Features(list(psms1 = psms, psms2 = psms), metadata = list())
  feat3 <- createOneToOneAssayLink(feat3, from = "psms1", to = "psms2", 
                                   fcol = "Var")
  al <- assayLink(feat3, "psms2")
  expect_identical(al@name, "psms2")
  expect_identical(al@from, "psms1")
  expect_identical(al@fcol, "Var")
  expect_identical(al@hits@from, 1:10)
  expect_identical(al@hits@to, 1:10)
  al@fcol <- "oneToOneID"
  expect_identical(assayLink(createOneToOneAssayLink(feat3, from = "psms1", 
                                                     to = "psms2"), "psms2"),
                   al)
  # Expect error: number of rows don't match
  expect_error(createOneToOneAssayLink(feat1, from = "psms", to = "peptides"))
  
})