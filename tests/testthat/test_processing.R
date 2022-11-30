data(feat1)

test_that("function: .applyTransformation", {
    data("feat2")
    featlog <- .applyTransformation(feat2, i = 1:2, 
                                    name = paste0("log", 1:2),
                                    .FUN = logTransform)
    ## Check the transformation is correctly applied
    expect_identical(featlog[["log1"]],
                     logTransform(featlog[[1]]))
    expect_identical(featlog[["log2"]],
                     logTransform(featlog[[2]]))
    ## Check the AssayLinks are correct
    expect_identical(featlog@assayLinks[["log1"]]@name, "log1")
    expect_identical(featlog@assayLinks[["log1"]]@from, "assay1")
    ## One-to-one mapping
    expect_identical(featlog@assayLinks[["log1"]]@hits@from,
                     featlog@assayLinks[["log1"]]@hits@to) 
    
    ## Test errors
    ## i and name have different lengths
    expect_error(.applyTransformation(feat1, i = 1, 
                                      name = letters[1:2],
                                      .FUN = logTransform),
                 regexp = "length.name. is not TRUE")
})

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

test_that("test transformation of multiple assays", {
    data("feat2")
    
    ## logTransform
    feat2_log <- logTransform(feat2, paste0("assay", 1:3), 
                              paste0("log", 1:3))
    for (i in 1:3) {
        expect_identical(feat2_log[[paste0("log", i)]],
                         logTransform(feat2_log[[paste0("assay", i)]]))
    }
    ## scaleTransform
    feat2_scale <- scaleTransform(feat2, paste0("assay", 1:3), 
                                  paste0("scale", 1:3))
    for (i in 1:3) {
        expect_identical(feat2_scale[[paste0("scale", i)]],
                         scaleTransform(feat2_scale[[paste0("assay", i)]]))
    }
    ## normalize
    feat2_norm <- normalize(feat2, paste0("assay", 1:3), 
                            paste0("norm", 1:3), method = "center.median")
    for (i in 1:3) {
        expect_identical(feat2_norm[[paste0("norm", i)]],
                         normalize(feat2_norm[[paste0("assay", i)]], 
                                   method = "center.median"))
    }
    ## sweep
    feat2_sweep <- sweep(feat2, i = paste0("assay", 1:3), 
                         name = paste0("sweep", 1:3), 
                         FUN = "-", MARGIN = 1, STATS = 1)
    for (i in 1:3) {
        expect_identical(feat2_sweep[[paste0("sweep", i)]],
                         sweep(feat2_sweep[[paste0("assay", i)]], 
                               FUN = "-", MARGIN = 1, STATS = 1))
    }
})
