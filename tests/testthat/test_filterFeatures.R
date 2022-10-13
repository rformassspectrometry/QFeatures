test_that("filterFeatures", {
    ## Prepare data
    data(feat1)    
    feat1 <- aggregateFeatures(feat1, 1, fcol = "Sequence", name = "peptides", fun = colMedians)
    feat1 <- aggregateFeatures(feat1, 2, fcol = "Protein", name = "proteins", fun = colMedians)
    
    ## Test character filters
    filter1 <- expect_message(filterFeatures(feat1, ~  location == "Mitochondrion"))
    filter2 <- expect_message(filterFeatures(feat1, ~startsWith(location, "Mito")))
    filter3 <- expect_message(filterFeatures(feat1, VariableFilter("location", "Mitochondrion")))
    filter4 <- expect_message(filterFeatures(feat1, VariableFilter("location", "unknown", condition = "!=")))
    filter5 <- expect_message(filterFeatures(feat1, VariableFilter("location", "unknown", condition = "==", not = TRUE)))
    filter6 <- expect_message(filterFeatures(feat1, ~ location != "unknown"))
    filter7 <- expect_message(filterFeatures(feat1, VariableFilter("location", "ochon", condition = "contains")))
    expect_equal(filter1, filter2)
    expect_equal(filter1, filter3)
    expect_equal(filter1, filter4)
    expect_equal(filter1, filter5)
    expect_equal(filter1, filter6)
    expect_equal(filter1, filter7)
    expect_identical(lengths(filter1), c(6L, 2L, 1L))
    
    ## Test filter stored in variable
    target <- "Mitochondrion"
    filter8 <- expect_message(filterFeatures(feat1, ~  location == target))
    filter9 <- expect_message(filterFeatures(feat1, VariableFilter("location", target)))
    expect_equal(filter1, filter8)
    expect_equal(filter8, filter9)
    ## Test filter stored in variable within function
    runfilter <- function() {
        target2 <- "Mitochondrion"
        expect_message(filterFeatures(feat1, ~  location == target2))
    }
    expect_equal(filter8, runfilter())
    
    ## Test numerical filters
    filter1 <- expect_message(filterFeatures(feat1, VariableFilter("pval", 0.03, "<=")))
    filter2 <- expect_message(filterFeatures(feat1, ~ pval <= 0.03))
    expect_equal(filter1, filter2)
    
    ## Test no match filters
    expect_true(all(dims(filterFeatures(feat1, VariableFilter("location", "not")))[1, ] == 0))
    expect_true(all(dims(filterFeatures(feat1, ~ location == "not"))[1, ] == 0))
    expect_true(all(dims(filterFeatures(feat1, ~ is.na(pval)))[1, ] == 0))
    
    ## Test fraud filters
    expect_error(VariableFilter("pval", TRUE, "<="))
    expect_error(VariableFilter("location", TRUE, "!="))
})  


test_that("filterFeatures with NAs", {
    ## Prepare data
    data(feat1)
    rowData(feat1[[1]])[1, "location"] <- NA
    flt <- VariableFilter("location", "Mitochondrion")
    
    ## Default is na.rm = FALSE
    res1 <-  expect_message(filterFeatures(feat1, flt))
    res2 <-  expect_message(filterFeatures(feat1, flt, na.rm = FALSE))
    expect_identical(res1, res2)
    expect_true(is.na(rowData(res1[[1]])[1, "location"]))
    
    ## Removing the row with the NA
    res3 <-  expect_message(filterFeatures(feat1, flt, na.rm = TRUE))
    expect_identical(nrow(assay(res1, 1)), nrow(assay(res3, 1)) + 1L)
    rd1 <- rowData(res1[[1]])
    rd3 <- rowData(res3[[1]])
    expect_identical(rd1[-1, ], rd3)

    ## VariableFilter and formula return identical results
    res4 <- expect_message(filterFeatures(feat1, ~ location == "Mitochondrion", na.rm = FALSE))
    res5 <- expect_message(filterFeatures(feat1, ~ location == "Mitochondrion", na.rm = TRUE))
    expect_identical(res2, res4)
    expect_identical(res3, res5)
})

test_that("filterFeatures with missing filter variables", {
    ## Prepare data
    data(feat1)    
    feat1 <- aggregateFeatures(feat1, 1, fcol = "Sequence", name = "peptides", fun = colMedians)
    feat1 <- aggregateFeatures(feat1, 2, fcol = "Protein", name = "proteins", fun = colMedians)
    
    ## Using Variable filter
    
    ## Test: 1 variable, missing in some assays, remove features for 
    ## missing variables
    expect_message(
        filt <- filterFeatures(feat1, VariableFilter("pval", 0.03, "<=")),
        regexp = "pval.*1 out of 3.*keep['] argument"
    )
    expect_identical(lengths(filt), c(3L, 0L, 0L))
    ## Test: 1 variable, missing in some assays, keep features for 
    ## missing variables
    expect_message(
        filt <- filterFeatures(feat1, VariableFilter("pval", 0.03, "<="), keep = TRUE),
        regexp = "pval.*1 out of 3.*keep['] argument"
    )
    expect_identical(lengths(filt), c(3L, 3L, 2L))
    ## Test: 1 variable, missing in all assays
    expect_error(filterFeatures(feat1, ~ foo),
                 regexp = "foo.*absent")
    ## Variable filter only allows filtering 1 variable at a time
    
    ## Using formula filter
    
    ## Test: 1 variable, missing in some assays, remove features for 
    ## missing variables
    expect_message(
        filt <- filterFeatures(feat1, ~ pval <= 0.03),
        regexp = "pval.*1 out of 3.*keep['] argument"
    )
    expect_identical(lengths(filt), c(3L, 0L, 0L))
    ## Test: 1 variable, missing in some assays, keep features for 
    ## missing variables
    expect_message(
        filt <- filterFeatures(feat1, ~ pval <= 0.03, keep = TRUE),
        regexp = "pval.*1 out of 3.*keep['] argument"
    )
    expect_identical(lengths(filt), c(3L, 3L, 2L))
    ## Test: 1 variable, missing in all assays
    expect_error(filterFeatures(feat1, ~ foo),
                 regexp = "foo.*absent")
    
    ## Test: 2 variables, 1 present in all assays and 1 in some, remove
    ## features for missing variables
    expect_message(
        filt <- filterFeatures(feat1, ~ pval <= 0.03 & grepl("Mito", location)),
        regexp = "pval.*1 out of 3.*location.*3 out of 3.*keep['] argument"
    )
    expect_identical(lengths(filt), c(2L, 0L, 0L))
    ## Test: 2 variables, 1 present in all assays and 1 in some, keep
    ## features for missing variables
    expect_message(
        filt <- filterFeatures(feat1, ~ pval <= 0.03 & grepl("Mito", location), keep = TRUE),
        regexp = "pval.*1 out of 3.*location.*3 out of 3.*keep['] argument"
    )
    expect_identical(lengths(filt), c(2L, 3L, 2L))
    ## Test: 2 variable, 1 missing in all assays
    expect_error(filterFeatures(feat1, ~ pval <= 0.03 & foo),
                 regexp = "foo.*absent")
    ## Test: 2 variable, 2 missing in all assays
    expect_error(filterFeatures(feat1, ~ foo & bar),
                 regexp = "foo['][,] [']bar.*absent")
})  
