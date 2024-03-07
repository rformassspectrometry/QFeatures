##' @title QFeatures from tabular data
##'
##' @description
##'
##' These functions convert tabular data into dedicate data
##' objets. The [readSummarizedExperiment()] function takes a
##' `data.frame` and converts it into a [SummarizedExperiment] object.
##' The [readQFeatures()] function takes a `data.frame` and converts
##' it into a [QFeatures] object. Two use-cases exist here:
##'
##' - The single-assay case will generate a [QFeatures] object with a
##'   single [SummarizedExperiment] assay containing all features of
##'   the input table.
##'
##' - The multi-assay case will generate a [QFeatures] object with
##'   multiple [SummarizedExperiment] assays, resulting from splitting
##'   the input table.
##'
##' @param assayData A `data.frame`, or any object that can be coerced
##'     to a `data.frame`, holding the quantitative assay.
##'
##' @param colAnnoation The type of this parameter will define whether
##'     the resulting [QFeaures] object will contain a single or
##'     multiple assays.
##'
##'     For the single-assay case, a `numeric` indicating the indices
##'     of the columns to be used as expression values, or a
##'     `character` indicating the names of the columns.
##'
##'     For the multi-assay case, a `data.frame` or any object that
##'     can be coerced to a `data.frame`. It is expected to contain
##'     all the sample meta information. Required fields are the
##'     acquisition batch (given by `batchCol`) and the acquisition
##'     channel within the batch (e.g. TMT channel, given by
##'     `channelCol`). Additional fields (e.g. sample type,
##'     acquisition date,...) are allowed and will be stored as sample
##'     meta data.
##'
##' @param fnames For the single-assay case, an optional
##'     `character(1)` or `numeric(1)` indicating the column to be
##'     used as feature names.
##'
##' @param name For the single-assay case, an optional `character(1)`
##'     to name the assay in the `QFeatures` object. If not set,
##'     `features` is used.
##'
##' @param batchCol For the multi-assay case, a `numeric(1)` or
##'     `character(1)` pointing to the column of `assayData` and
##'     `colAnnotation` that contain the batch names. Make sure that
##'     the column name in both table are either identical and
##'     syntactically valid (if you supply a `character`) or have the
##'     same index (if you supply a `numeric`). Note that characters
##'     can be converted to syntactically valid names using
##'     `make.names`
##'
##' @param channelCol For the multi-assay case, a `numeric(1)` or
##'     `character(1)` pointing to the column of `colData` that
##'     contains the column names of the quantitative data in
##'     `featureData` (see Example).
##'
##' @param suffix For the multi-assay case, a `character()` giving the
##'     suffix of the column names in each assay. Sample/single-cell
##'     (column) names are automatically generated using: batch name +
##'     sep + suffix. Make sure suffix contains unique character
##'     elements. The length of the vector should equal the number of
##'     quantification channels.  If `NULL` (default), the suffix is
##'     derived from the the names of the quantification columns in
##'     `assayData`.
##'
##' @param sep A `character(1)` that is inserted between the assay
##'     name and the `suffix` (see `suffix` argument for more
##'     details).
##'
##' @param removeEmptyCols A `logical(1)`. If true, the function will
##'     remove in each batch the columns that contain only missing
##'     values.
##'
##' @param verbose A `logical(1)` indicating whether the progress of
##'     the data reading and formatting should be printed to the
##'     console. Default is `TRUE`.
##'
##' @return An instance of class [QFeatures] or
##'     [SummarizedExperiment]. For the former, the expression data of
##'     each batch is stored in a separate assay as a
##'     [SummarizedExperiment] object.
##'
##' @author Laurent Gatto, Christophe Vanderaa
##'
##' @importFrom methods new validObject
##' @import SummarizedExperiment
##'
##' @seealso The [QFeatures] class for an example on how to use
##'     `readQFeatures` and how to further manipulate the resulting
##'     data.
##'
##' @md
##'
##' @name readQFeatures
##' @aliases readSummarizedExperiment
##' @aliases readQFeatures
##' @aliases readQFeatures,data.frame,data.frame
##' @aliases readQFeatures,data.frame,vector
##' @export
##'
##' @examples
##'
##' ###################################
##' ## Single-assay case.
##'
##' ## Load a data.frame with PSM-level data
##' data(hlpsms)
##'
##' ## Create the QFeatures object
##' fts2 <- readQFeatures(hlpsms, colAnnotation = 1:10, name = "psms")
##' fts2
##'
##' ###################################
##' ## Multi-assay case.
##' ## See scp::readSCP()
NULL

## Simple case, with a single table to populate one QFeatures
## assay.

##' @export
setMethod("readQFeatures", c("data.frame", "vector"),
          function(assayData, colAnnotation,
                   fnames, name = NULL)
              .readQFeatures1(assayData, colAnnotation,
                              fnames, name = NULL))

## Second case, with a single table to populate multiple QFeatures
## assay. Only from a data.frame (not a file name), handled by
## .readQFeatures2(). The second argument is the colData.

##' @export
setMethod("readQFeatures", c("data.frame", "data.frame"),
          function(assayData, colAnnotation,
                   batchCol, channelCol, suffix = NULL, sep = "",
                   removeEmptyCols = FALSE, verbose = TRUE)
              .readQFeatures2(featureData, colData, batchCol, channelCol,
                              suffix = NULL, sep = "", removeEmptyCols = FALSE,
                              verbose = TRUE))
##' @export
##'
##' @rdname readQFeatures
##'
##' @importFrom utils read.csv
##'
##' @param ... Further arguments that can be passed on to [read.csv()]
##'     except `stringsAsFactors`, which is always `FALSE`.
readSummarizedExperiment <- function(assayData, colAnnotation,
                                     fnames, ...) {
    ecol <- colAnnotation ## still use former arg name
    if (is.data.frame(assayData)) xx <- assayData
    else {
        args <- list(...)
        args$file <- assayData
        if ("rownames" %in% names(args)) {
            if (missing(fnames)) fnames <- args$rownames
            args$rownames <- NULL
        }
        args$stringsAsFactors <- FALSE
        xx <- do.call(read.csv, args)
    }
    if (is.character(ecol) || is.factor(ecol)) {
        ecol0 <- ecol
        ecol <- match(ecol0, colnames(xx))
        if (any(is.na(ecol)))
            stop("Column identifiers ",
                 paste(ecol0[is.na(ecol)], collapse = ", "),
                 " not recognised among\n",
                 paste(colnames(xx), paste = ", "))
    } else if (is.logical(ecol)) {
        if (length(ecol) != length(xx))
            stop("Length of 'colAnnotation' and 'assayData' do not match.")
        ecol <- which(ecol)
    }
    assay <- as.matrix(xx[, ecol, drop = FALSE])
    fdata <- DataFrame(xx[, -ecol, drop = FALSE])

    if (!missing(fnames)) {
        fnames <- fnames[1]
        if (is.numeric(fnames))
            fnames <- colnames(xx)[fnames]
        if (is.na(match(fnames, colnames(xx))))
            stop(fnames, " not found among\n",
                 paste(colnames(xx), paste = ", "))
        rownames(fdata) <- rownames(assay) <- fdata[, fnames]
    } else {
        rownames(fdata) <- rownames(assay) <- seq_len(nrow(assay))
    }
    SummarizedExperiment(assay, rowData = fdata)
}


##' @export
##'
##' @rdname QFeatures-class
##'
##' @param assayLinks An optional [AssayLinks] object.
QFeatures <- function(..., assayLinks = NULL) {
    ans <- MultiAssayExperiment(...)
    if (isEmpty(ans)) assayLinks <- AssayLinks()
    else {
        if (is.null(assayLinks))
            assayLinks <- AssayLinks(names = names(ans))
    }
    new("QFeatures",
        ExperimentList = ans@ExperimentList,
        colData = ans@colData,
        sampleMap = ans@sampleMap,
        metadata = ans@metadata,
        assayLinks = assayLinks)
}


.readQFeatures1 <- function(assayData, colAnnotation,
                            fnames, name = NULL)  {
    se <- readSummarizedExperiment(assayData, colAnnotation,
                                   fnames)
    if (anyDuplicated(rownames(se))) {
        message("Making assay rownames unique.")
        rownames(se) <- make.unique(rownames(se))
    }
    cd <- DataFrame(row.names = colnames(se))
    if (is.null(name))
        name <- "features"
    el <- structure(list(se), .Names = name[1])
    al <- AssayLinks(AssayLink(name = name[1]))
    ans <- MultiAssayExperiment(el, colData = cd)
    new("QFeatures",
        ExperimentList = ans@ExperimentList,
        colData = ans@colData,
        sampleMap = ans@sampleMap,
        metadata = ans@metadata,
        assayLinks = al)
}

.readQFeatures2 <- function(assayData, colAnnotation,
                            batchCol, channelCol, suffix = NULL,
                            sep = "", removeEmptyCols = FALSE,
                            verbose = TRUE) {
    ## Check the batch column name
    if (!identical(make.names(batchCol), batchCol))
        stop("'batchCol' is not a syntactically valid column name. ",
             "See '?make.names' for converting the column names to ",
             "valid names, e.g. '", batchCol, "' -> '",
             make.names(batchCol), "'")
    colData <- as.data.frame(colData)
    ## Get the column contain the expression data
    ecol <- unique(colData[, channelCol])
    ## Get the sample suffix
    if (is.null(suffix))
        suffix <- ecol
    ## Create the SummarizedExperiment object
    if (verbose) message("Loading data as a 'SummarizedExperiment' object")
    se <- readSummarizedExperiment(assayData,
                                   colAnnotation)
    if (is.null(list(...)$row.names))
        rownames(se) <- paste0("PSM", seq_len(nrow(se)))
    ## Check the link between colData and se
    mis <- !rowData(se)[, batchCol] %in% colData[, batchCol]
    if (any(mis)) {
        warning("Missing metadata. The features are removed for ",
                paste0(unique(rowData(se)[mis, batchCol]), collapse = ", "))
        se <- se[!mis, ]
    }
    ## Split the SingleCellExperiment object by batch column
    if (verbose) message("Splitting data based on '", batchCol, "'")
    se <- .splitSE(se, f = batchCol)
    ## Clean each element in the data list
    for (i in seq_along(se)) {
        ## Add unique sample identifiers
        colnames(se[[i]]) <- paste0(names(se)[[i]], sep, suffix)
        ## Remove the columns that are all NA
        if (removeEmptyCols) {
            sel <- colSums(is.na(assay(se[[i]]))) != nrow(se[[i]])
            se[[i]] <- se[[i]][, sel]
        }
    }
    if (verbose) message("Formatting sample metadata (colData)")
    ## Create the colData
    cd <- DataFrame(row.names = unlist(lapply(se, colnames)))
    rownames(colData) <- paste0(colData[, batchCol], sep, suffix)
    cd <- cbind(cd, colData[rownames(cd), ])
    ## Store the data as a QFeatures object
    if (verbose) message("Formatting data as a 'QFeatures' object")
    QFeatures(experiments = se, colData = cd)
}
