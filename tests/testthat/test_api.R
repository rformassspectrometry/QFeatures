## version 1.17.2
api <- c("addAssay", "addAssayLink", "addAssayLinkOneToOne", "adjacencyMatrix",
         "adjacencyMatrix<-", "aggcounts", "aggregateFeatures", "assayLink",
         "AssayLink", "assayLinks", "AssayLinks", "coerce",
         "countUniqueFeatures", "dims", "display", "dropEmptyAssays",
         "expandDataFrame", "filterFeatures", "filterNA", "impute", "infIsNA",
         "joinAssays", "logTransform", "longFormat", "longForm", "ncols", "nNA",
         "normalize", "nrows", "QFeatures", "rbindRowData", "readQFeatures",
         "readQFeaturesFromDIANN", "readSummarizedExperiment",
         "reduceDataFrame", "removeAssay", "replaceAssay", "rowData<-",
         "rowDataNames", "scaleTransform", "selectRowData", "show",
         "subsetByFeature", "sweep", "unfoldDataFrame", "updateObject",
         "VariableFilter", "zeroIsNA")

test_that("API hasn't changed", {
    current_api <- sort(ls(pos = "package:QFeatures"))
    expect_identical(length(current_api), length(api))
    expect_identical(current_api, sort(api))
})
