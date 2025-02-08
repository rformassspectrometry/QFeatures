data(feat3)
qf <- feat3[, , 1:3]

test_that("createPrecursorId work", {
    expect_error(createPrecursorId(feat3, name = "Sequence"),
                 "'Sequence' already exists.")
    expect_error(createPrecursorId(feat3,
                                   fcols = c("Sequence", "Var")),
                 "'fcols' not found in some assays.")
    res <- createPrecursorId(qf,
                             fcols = c("Sequence", "Var"))
    expect_identical(lengths(rowDataNames(qf)) + 1L,
                     lengths(rowDataNames(res)))
    expect_true(all(sapply(rowDataNames(res), "[[", 6) == "Precursor.Id"))
})
