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
    al <- assayLink(feat1_log, "log_psms")
    expect_identical(al, 
                     .createAssayLinkOneToOne(seFrom = feat1_log[["psms"]],
                                              seTo = feat1_log[["log_psms"]],
                                              nameFrom = "psms",
                                              nameTo = "log_psms"))
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
    al <- assayLink(feat1_scaled, "scaled_psms")
    expect_identical(al, 
                     .createAssayLinkOneToOne(seFrom = feat1_scaled[["psms"]],
                                              seTo = feat1_scaled[["scaled_psms"]],
                                              nameFrom = "psms",
                                              nameTo = "scaled_psms"))
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
    al <- assayLink(feat1_norm, "norm_psms")
    expect_identical(al, 
                     .createAssayLinkOneToOne(seFrom = feat1_norm[["psms"]],
                                              seTo = feat1_norm[["norm_psms"]],
                                              nameFrom = "psms",
                                              nameTo = "norm_psms"))
})

test_that("function: all normalize methods", {
    data(hlpsms)
    fts <- readFeatures(hlpsms[1:5000,], ecol = 1:10, name = "psms")
    se <- fts[[1]]
    for (.method in MsCoreUtils::normalizeMethods()) {
        se_norm <- normalize(se, method = .method)
        feat_norm <- normalize(fts, 1, method = .method)
        expect_identical(se_norm, feat_norm[["normAssay"]])
    }
})
