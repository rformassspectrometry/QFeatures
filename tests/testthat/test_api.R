## version 1.19.1
api <- c("addAssay", "addAssayLink", "addAssayLinkOneToOne",
         "adjacencyMatrix", "adjacencyMatrix<-", "aggcounts",
         "aggregateFeatures", "assayLink", "AssayLink", "assayLinks",
         "AssayLinks", "coerce", "countUniqueFeatures", "dims",
         "display", "dropEmptyAssays", "expandDataFrame",
         "filterFeatures", "filterNA", "getQFeaturesType", "impute", "infIsNA",
         "joinAssays", "logTransform", "longFormat", "longForm", "ncols", "nNA",
         "normalize", "nrows", "QFeatures", "rbindRowData",
         "readQFeatures", "readQFeaturesFromDIANN",
         "readSummarizedExperiment", "reduceDataFrame", "removeAssay",
         "replaceAssay", "rowData<-", "rowDataNames",
         "scaleTransform", "selectRowData", "setQFeaturesType", "show",
         "subsetByFeature", "sweep", "unfoldDataFrame", "updateObject",
         "validQFeaturesTypes", "VariableFilter", "zeroIsNA")

test_that("API hasn't changed", {
    current_api <- sort(ls(pos = "package:QFeatures"))
    expect_identical(length(current_api), length(api))
    expect_identical(current_api, sort(api))
})