##' @export
combineFeatures <- function(object,
                            i,
                            fcol,
                            name = NULL,
                            fun = median,
                            ... 
                            ) {
    stopifnot(inherits(object, "Features"))
    if (isEmpty(object))
        return(object)
    if (missing(i))
        i <- main_assay(object)
    .assay <- assay(object, i)
    .featureData <- featureData(object, i)
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

    .assay <- combine_assay(.assay, groupBy, fun, ...)
    .featureData <- reduce_DataFrame(.featureData, .featureData[[fcol]],
                                     simplify = TRUE, count = TRUE)
    
    new_fs <- new("FeatureSet",
                  assay = .assay,
                  featureData = .featureData[rownames(.assay), ],
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


combine_assay <- function(assay, groupBy, fun, ...) 
    do.call(rbind,
            by(assay, groupBy,
               function(.x) apply(.x, 2, fun, ...)))

