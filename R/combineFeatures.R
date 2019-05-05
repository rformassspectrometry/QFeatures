combineFeatures <- function(object,
                            i,
                            fcol,
                            name = "newAssay",
                            fun = median,
                            ... 
                            ) {
    if (isEmpty(object))
        return(object)
    if (missing(i))
        i <- main_assay(object)
    .assay <- assay(object, i)
    .rowdata <- rowData(object[[i]])
    if (missing(fcol))
        stop("Require either 'groupBy' or 'fcol'.")
    stopifnot(fcol %in% names(.rowdata))
    groupBy <- .rowdata[[fcol]]

    if (anyNA(.assay)) {
        msg <- paste("Your data contains missing values.",
                     "Please read the relevant section in the",
                     "combineFeatures manual page for details the",
                     "effects of missing values on data aggregation.")
        message(paste(strwrap(msg), collapse = "\n"))
    }

    .assay <- combine_assay(.assay, groupBy, fun, ...)
    .featureData <- Features::reduce_DataFrame(.rowdata, .rowdata[[fcol]],
                                     simplify = TRUE, count = TRUE)
    
    se <- SummarizedExperiment(.assay,
                               rowData = .featureData[rownames(.assay), ])
    assayLinks <- AssayLinks(name = name,
                             from = ifelse(is.character(i), i, names(object)[i]),
                             fcol = fcol)
    object <- addAssay(object, se, name = name, assayLinks = assayLinks)

    if (validObject(object))
        object
}


combine_assay <- function(assay, groupBy, fun, ...) 
    do.call(rbind,
            by(assay, groupBy,
               function(.x) apply(.x, 2, fun, ...)))
