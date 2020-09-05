data(feat1)

test_that("function: logTransform", {
    se <- feat1[[1]]
    se_log <- logTransform(se)
    feat1_log <- logTransform(feat1, 1)
    expect_identical(names(feat1_log), c("psms", "logAssay"))
    feat1_log <- logTransform(feat1, 1, "log_psms")
    expect_identical(names(feat1_log), c("psms", "log_psms"))
    expect_identical(se_log, feat1_log[[2]])
    e <- log2(assay(se))
    expect_identical(assay(se_log), e)
    expect_identical(assay(logTransform(se, base = 10)),
                     log10(assay(se)))
    ## Check the one-to-one link
    feat_log_sub <- feat1_log["PSM1", , ]
    expect_identical(log2(assay(feat_log_sub, 1)),
                     assay(feat_log_sub, 2))
})


test_that("function: scaleTransform", {
    se <- feat1[[1]]
    se_scaled <- scaleTransform(se)
    feat1_scaled <- scaleTransform(feat1, 1)
    expect_identical(names(feat1_scaled), c("psms", "scaledAssay"))
    feat1_scaled <- scaleTransform(feat1, 1, "scaled_psms")
    expect_identical(names(feat1_scaled), c("psms", "scaled_psms"))
    expect_identical(se_scaled, feat1_scaled[[2]])
    e <- scale(assay(se), center = TRUE, scale = TRUE)
    attr(e, "scaled:center") <-
        attr(e, "scaled:scale") <- NULL
    expect_identical(assay(se_scaled), e)
    ## Check the one-to-one link
    feat_scaled_sub <- feat1_scaled["PSM1", , ]
    expect_identical(dimnames(assay(feat_scaled_sub, 1)),
                     dimnames(assay(feat_scaled_sub, 2)))
})


test_that("function: normalize", {
    se <- feat1[[1]]
    se_norm <- normalize(se, method = "max")
    feat1_norm <- normalize(feat1, 1, method = "max")
    expect_identical(names(feat1_norm), c("psms", "normAssay"))
    feat1_norm <- normalize(feat1, 1, "norm_psms", method = "max")
    expect_identical(names(feat1_norm), c("psms", "norm_psms"))
    expect_identical(se_norm, feat1_norm[[2]])
    e <- assay(se) / rowMax(assay(se))
    expect_identical(assay(se_norm), e)
    ## Check the one-to-one link
    feat_norm_sub <- feat1_norm["PSM1", , ]
    expect_identical(dimnames(assay(feat_norm_sub, 1)),
                     dimnames(assay(feat_norm_sub, 2)))
})

test_that("function: all normalize methods", {
    data(hlpsms)
    fts <- readQFeatures(hlpsms, ecol = 1:10, name = "psms")
    se <- fts[[1]]
    for (.method in MsCoreUtils::normalizeMethods()) {
        se_norm <- normalize(se, method = .method)
        feat_norm <- normalize(fts, 1, method = .method)
        expect_identical(se_norm, feat_norm[["normAssay"]])
    }
})



test_that("function: sweep", {
    se <- feat1[[1]]
    cmeds <- colMedians(assay(se))
    e <- sweep(assay(se), 2, cmeds)
    sse <- sweep(se, 2, cmeds)
    sse2 <- QFeatures:::sweepSE(se, 2, cmeds)
    expect_identical(e, assay(sse))
    expect_identical(e, assay(sse2))
    sfeat1 <- sweep(feat1, MARGIN = 2, STATS = cmeds, i = 1)
    expect_identical(names(sfeat1), c("psms", "sweptAssay"))
    expect_identical(e, assay(sfeat1[["sweptAssay"]]))
})
