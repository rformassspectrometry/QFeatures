combineFeatures <- function(object,
                            i,
                            fcol,
                            name = NULL,
                            method = c("mean",
                                       "median",
                                       "weighted.mean",
                                       "sum",
                                       "medpolish",
                                       "robust",
                                       "iPQF",
                                       "NTR"),
                            redundancy.handler = c("unique", "multiple"),
                            cv = FALSE,
                            cv.norm = "sum",
                            ... ## further arguments to method
                            ) {
    stopifnot(inherits(object, "Features"))
    if (isEmpty(object))
        return(object)
    if (missing(i))
        i <- main_assay(object)
    .assay <- assay(object, i)
    .featureData <- featureData(object, i)
    if (is.character(method))
        method <- match.arg(method)
    if (missing(fcol))
        stop("Require either 'groupBy' or 'fcol'.")
    stopifnot(fcol %in% names(.featureData))
    groupBy <- .featureData[[fcol]]

    if (anyNA(.assay)) {
        msg <- paste("Your data contains missing values.",
                     "Please read the relevant section in the",
                     "combineFeatures manual page for details the",
                     "effects of missing values on data aggregation.")
        message(paste(strwrap(msg), collapse = "\n"))
    }

    ## .assay <- combine_features_assay(.assay, groupBy, ...)
    ## .featureData <- combine_features_featureData(.featureData, groupBy, ...)
    new_fs <- new("FeatureSet",
                  assay = .assay,
                  featureData = .featureData,
                  id = get_next_featureSet_id(object),
                  from = slot(object[[i]], "id"))
    .features <- features(object)
    .features <- append(.features, new_fs)
    if (!is.null(name))
        names(.features)[length(.features)] <- name
    object@listData <- .features
    if (validObject(object))
        object
}
