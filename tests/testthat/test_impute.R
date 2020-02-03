data(se_na2)
x <- assay(se_na2)
randna  <-  rowData(se_na2)$randna

test_that("impute: mandatory method", {
    data(se_na2)
    expect_error(impute(se_na2))
    expect_error(impute(se_na2, method = "not"))
})

test_that("impute: absence of missing values", {
    data(se_na2)
    se_imp <- impute(se_na2, method = "knn")
    se_imp_2 <- impute(se_imp, method = "knn")
    expect_identical(se_imp, se_imp_2)    
})

test_that("impute,SummarizedExperiment", {
    data(se_na2)
    x <- assay(se_na2)
    se_imp <- impute(se_na2, method = "knn")
    x_imp <- MsCoreUtils::impute_matrix(x, method = "knn")
    expect_identical(x_imp, assay(se_imp))    
})

test_that("impute,Features", {
    data(se_na2)
    ft <- Features(list(se_na2 = se_na2), colData = colData(se_na2))
    x <- assay(se_na2)
    ft_imp <- impute(ft, method = "MinDet")
    x_imp <- MsCoreUtils::impute_matrix(x, method = "MinDet")
    expect_identical(x_imp, assay(ft_imp[[1]]))
})
