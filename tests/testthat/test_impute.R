data(se_na2)
x <- assay(se_na2)
randna  <-  rowData(se_na2)$randna
ft <- QFeatures(list(se_na2 = se_na2), colData = colData(se_na2))

test_that("impute: mandatory method", {
    data(se_na2)
    expect_error(impute(se_na2))
    expect_error(impute(se_na2, method = "not"))
})

test_that("impute: test warning", {
    data(se_na2)
    ## Get the number of rows with more than 50% missing
    nnarows <- sum(rowMeans(is.na(assay(se_na2))) > 0.5)
    ## Test warning
    expect_warning(impute(se_na2, method = "knn"),
                   regexp = paste0(nnarows, " rows with more than 50 %"))
})

test_that("impute: absence of missing values", {
    data(se_na2)
    se_imp <- expect_warning(impute(se_na2, method = "knn"))
    se_imp_2 <- impute(se_imp, method = "knn")
    expect_identical(se_imp, se_imp_2)
})

test_that("impute,SummarizedExperiment", {
    data(se_na2)
    x <- assay(se_na2)
    se_imp <- expect_warning(impute(se_na2, method = "knn"))
    x_imp <- expect_warning(MsCoreUtils::impute_matrix(x, method = "knn"))
    expect_identical(x_imp, assay(se_imp))
})

test_that("impute,QFeatures: test errors", {
    ## i is mandatory
    expect_error(impute(ft, method = "zero"), 
                 regexp = "Provide index")
    ## Only one sample can be imputed at once
    expect_error(impute(ft, i = 1:2, method = "zero"), 
                 regexp = "one assay")
})

test_that("impute,QFeatures", {
    ## Check imputation
    ft_imp <- impute(ft, i = "se_na2", method = "MinDet")
    x_imp <- MsCoreUtils::impute_matrix(x, method = "MinDet")
    ## Check the new assay was correctly imputed
    expect_identical(x_imp, assay(ft_imp[["imputedAssay"]]))
    ## Check the assay was added
    expect_true(length(ft_imp) == (length(ft) + 1))
})
