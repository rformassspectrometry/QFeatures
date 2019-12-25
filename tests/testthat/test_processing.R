data(feat1)

test_that("function: logTransform", {
    se <- feat1[[1]]
    se_log <- logTransform(se)
    feat1_log <- logTransform(feat1)
    expect_identical(se_log, feat1_log[[1]])
    e <- log2(assay(se))
    expect_identical(assay(se_log), e)
    expect_identical(assay(logTransform(se, base = 10)),
                     log10(assay(se)))
})


test_that("function: scaleTransform", {
    se <- feat1[[1]]
    se_scaled <- scaleTransform(se)
    feat1_scaled <- scaleTransform(feat1)
    expect_identical(se_scaled, feat1_scaled[[1]])
    e <- scale(assay(se), center = TRUE, scale = TRUE)
    attr(e, "scaled:center") <-
        attr(e, "scaled:scale") <- NULL
    expect_identical(assay(se_scaled), e)    
})


test_that("function: normalize", {
    se <- feat1[[1]]
    se_norm <- normalize(se, method = "max")
    feat1_norm <- normalize(feat1, method = "max")
    expect_identical(se_norm, feat1_norm[[1]])
    e <- assay(se) / rowMax(assay(se))
    expect_identical(assay(se_norm), e)    
})

test_that("function: all normalize methods", {
    data(hlpsms)
    fts <- readFeatures(hlpsms[1:5000,], ecol = 1:10, name = "psms")
    se <- fts[[1]]
    for (.method in normalizeMethods()) {
        se_norm <- normalize(se, method = .method)
        feat_norm <- normalize(fts, method = .method)
        expect_identical(se_norm, feat_norm[[1]])
    }
})
