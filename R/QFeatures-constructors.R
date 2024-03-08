#' @title QFeatures from tabular data
##'
##' @description
##'
##' These functions convert tabular data into dedicate data
##' objets. The [readSummarizedExperiment()] function takes a
##' `data.frame` and converts it into a [SummarizedExperiment] object.
##' The [readQFeatures()] function takes a `data.frame` and converts
##' it into a `QFeatures` object (see [QFeatures()] for details). Two
##' use-cases exist here:
##'
##' - The single-set case will generate a `QFeatures` object with a
##'   single [SummarizedExperiment] set containing all features of the
##'   input table.
##'
##' - The multi-set case will generate a `QFeatures` object with
##'   multiple [SummarizedExperiment] sets, resulting from splitting
##'   the input table.
##'
##' @param assayData A `data.frame`, or any object that can be coerced
##'     to a `data.frame`, holding the quantitative assay. For
##'     `readSummarizedExperiment()`, this can also be a
##'     `character(1)` to a filename.
##'
##' @param colAnnotation The type of this parameter will define
##'     whether the resulting `QFeaures` object will contain a single
##'     or multiple sets.
##'
##'     For the single-set case, a `numeric` indicating the indices of
##'     the columns to be used as expression values, or a `character`
##'     indicating the names of the columns.
##'
##'     For the multi-set case, a `data.frame` or any object that can
##'     be coerced to a `data.frame`. It is expected to contain all
##'     the sample meta information. Required fields are the
##'     acquisition batch (given by `batchCol`) and the acquisition
##'     channel within the batch (e.g. TMT channel, given by
##'     `channelCol`). Additional fields (e.g. sample type,
##'     acquisition date,...) are allowed and will be stored as sample
##'     meta data.
##'
##' @param fnames For the single- and multi-set cases, an optional
##'     `character(1)` or `numeric(1)` indicating the column to be
##'     used as feature names. If missing, in the latter case, the
##'     rows are named PSM1, PSM2, PSM3, ... Note that rownames must
##'     be unique in `QFeatures` sets.
##'
##' @param name For the single-set case, an optional `character(1)` to
##'     name the set in the `QFeatures` object. If not set, `psms`
##'     is used.
##'
##' @param batchCol For the multi-set case, a `numeric(1)` or
##'     `character(1)` pointing to the column of `assayData` and
##'     `colAnnotation` that contain the batch names. Make sure that
##'     the column name in both table are either identical and
##'     syntactically valid (if you supply a `character`) or have the
##'     same index (if you supply a `numeric`). Note that characters
##'     can be converted to syntactically valid names using
##'     `make.names`
##'
##' @param channelCol For the multi-set case, a `numeric(1)` or
##'     `character(1)` pointing to the column of `colData` that
##'     contains the column names of the quantitative data in
##'     `assayData` (see Example).
##'
##' @param suffix For the multi-set case, a `character()` giving the
##'     suffix of the column names in each set. Sample/single-cell
##'     (column) names are automatically generated using: batch name +
##'     sep + suffix. Make sure suffix contains unique character
##'     elements. The length of the vector should equal the number of
##'     quantification channels.  If `NULL` (default), the suffix is
##'     derived from the the names of the quantification columns in
##'     `assayData`.
##'
##' @param sep A `character(1)` that is inserted between the set
##'     name and the `suffix` (see `suffix` argument for more
##'     details).
##'
##' @param removeEmptyCols A `logical(1)`. If true, the function will
##'     remove in each batch the columns that contain only missing
##'     values.
##'
##' @param ecol Same as `colAnnotation` for the single-set
##'     case. Available for backwards compatibility. Default is
##'     `NULL`. If both `ecol` and `colAnnotation` are set, an error
##'     is thrown.
##'
##' @param verbose A `logical(1)` indicating whether the progress of
##'     the data reading and formatting should be printed to the
##'     console. Default is `TRUE`.
##'
##' @return An instance of class `QFeatures` or
##'     [SummarizedExperiment::SummarizedExperiment()]. For the
##'     former, the expression data of each batch is stored in a
##'     separate set as a
##'     [SummarizedExperiment::SummarizedExperiment()] object.
##'
##' @author Laurent Gatto, Christophe Vanderaa
##'
##' @importFrom methods new validObject
##' @import SummarizedExperiment
##'
##' @seealso
##'
##' - The `QFeatures` (see [QFeatures()]) class to read about how to
##'   manipulate the resulting `QFeatures` object.
##'
##' - The [readQFeaturesfromDIANN()] function to import DIA-NN
##'   quantitative data.
##'
##' @md
##'
##' @name readQFeatures
##' @aliases readSummarizedExperiment
##' @aliases readQFeatures
##' @aliases readQFeatures,data.frame,data.frame
##' @aliases readQFeatures,data.frame,vector
##' @aliases readQFeatures,missing,vector
##' @export
##'
##' @examples
##'
##' ######################################
##' ## Single-set case.
##'
##' ## Load a data.frame with PSM-level data
##' data(hlpsms)
##' hlpsms[1:10, c(1, 2, 10:11, 14, 17)]
##'
##' ## Create the QFeatures object
##' qf1 <- readQFeatures(hlpsms, colAnnotation = 1:10, name = "psms")
##' qf1
##' colData(qf1)
##'
##' ######################################
##' ## Single-set case using a data.frame
##'
##' ## All PSMs were acquired in the same acquisition
##' hlpsms$file <- "File1"
##' hlpsms[1:10, c(1, 2, 10:11, 14, 17, 29)]
##' (colann <- data.frame(file = rep("File1", 10),
##'                       Channel = names(hlpsms)[1:10]))
##'
##' qf2 <- readQFeatures(hlpsms, colAnnotation = colann,
##'                      batchCol = "file", channelCol = "Channel")
##' qf2
##' colData(qf2)
##'
##' ######################################
##' ## Multi-set case.
##'
##' ## Let's simulate 3 different files/batches for that same input
##' ## data.frame, and define a colAnnotation data.frame.
##'
##' hlpsms$file <- paste0("File", sample(1:3, nrow(hlpsms), replace = TRUE))
##' hlpsms[1:10, c(1, 2, 10:11, 14, 17, 29)]
##' (colann <- data.frame(file = rep(paste0("File", 1:3), each = 10),
##'                       Channel = rep(names(hlpsms)[1:10], 3)))
##'
##' qf3 <- readQFeatures(hlpsms, colAnnotation = colann,
##'                      batchCol = "file", channelCol = "Channel")
##' qf3
##' colData(qf3)
NULL

## Simple case, with a single table to populate one QFeatures set

##' @rdname readQFeatures
##' @export
setMethod("readQFeatures", c("data.frame", "missing"),
          function(assayData, colAnnotation,
                   fnames, name = NULL, ecol = NULL) {
              if (is.null(ecol))
                  stop("Please provide a 'colAnnotaion'.")
              stopifnot(is.vector(ecol))
              colAnnotation <- ecol
              .readQFeatures1(assayData, colAnnotation,
                              fnames, name)
          })

##' @rdname readQFeatures
##' @export
setMethod("readQFeatures", c("data.frame", "vector"),
          function(assayData, colAnnotation,
                   fnames, name = NULL)
              .readQFeatures1(assayData, colAnnotation,
                              fnames, name))



## Second case, with a single table to populate multiple QFeatures
## sets. Only from a data.frame, handled by .readQFeatures2().

##' @rdname readQFeatures
##' @export
setMethod("readQFeatures", c("data.frame", "data.frame"),
          function(assayData, colAnnotation,
                   batchCol, channelCol, suffix = NULL, sep = "",
                   removeEmptyCols = FALSE, fnames, verbose = TRUE,
                   ecol = NULL) {
              if (!is.null(ecol))
                  warning("Using 'colAnnotation' and ignoring 'ecol'.")
              .readQFeatures2(assayData, colAnnotation, batchCol,
                              channelCol, suffix, sep,
                              removeEmptyCols, fnames, verbose)
          })
##' @export
##'
##' @rdname readQFeatures
##'
##' @importFrom utils read.csv
##'
##' @param ... Further arguments that can be passed on to [read.csv()]
##'     except `stringsAsFactors`, which is always `FALSE`. Only
##'     applicable to `readSummarizedExperiment()`.
readSummarizedExperiment <- function(assayData, colAnnotation,
                                     fnames, ecol = NULL, ...) {
    if (missing(colAnnotation)) {
        if (is.null(ecol))
            stop("Please provide a 'colAnnotaion'.")
    } else {
        ## There is a colAnnotation
        if (!is.vector(colAnnotation))
            stop("'colAnnotation', must be a vector.")
        ecol <- colAnnotation
    }
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
        name <- "psms"
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
                            fnames, verbose = TRUE) {
    if (missing(batchCol))
        stop("Please provide a 'batchCol'.")
    if (missing(channelCol))
        stop("Please provide a 'channelCol'.")
    ## Check the batch column name
    if (!identical(make.names(batchCol), batchCol))
        stop("'batchCol' is not a syntactically valid column name. ",
             "See '?make.names' for converting the column names to ",
             "valid names, e.g. '", batchCol, "' -> '",
             make.names(batchCol), "'")
    colData <- as.data.frame(colAnnotation)
    ## Get the column contain the expression data
    ecol <- unique(colData[, channelCol])
    ## Get the sample suffix
    if (is.null(suffix))
        suffix <- ecol
    ## Create the SummarizedExperiment object
    if (verbose) message("Loading data as a 'SummarizedExperiment' object")
    se <- readSummarizedExperiment(assayData, ecol)
    if (missing(fnames)) {
        rownames(se) <- paste0("PSM", seq_len(nrow(se)))
    } else {
        stopifnot(fnames %in% names(rowData(se)))
        rownames(se) <- rowData(se)[, fnames]
    }
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



##' @title  Read DIA-NN output as a QFeatures objects
##'
##' @description
##'
##' This function takes the output tables from DIA-NN and converts them
##' into a `QFeatures` object.
##'
##' @param colData A `data.frame` or any object that can be coerced to
##'     a `data.frame`. `colData` is expected to contains all the
##'     sample annotations. We require the table to contain a column
##'     called `File.Name` that links to the `File.Name` in the DIA-NN
##'     report table. If `multiplexing = "mTRAQ"`, we require a second
##'     column called `Label` that links the label to the sample (the
##'     labels identified by DIA-NN can be retrieved from `Modified
##'     Sequence` column in the report table).
##'
##' @param reportData A `data.frame` or any object that can be coerced
##'     to a data.frame that contains the data from the `Report.tsv`
##'     file generated by DIA-NN.
##'
##' @param extractedData A data.frame or any object that can be coerced
##'     to a data.frame that contains the data from the `*_ms1_extracted.tsv`
##'     file generated by DIA-NN. This argument is optional and is
##'     only applicable for mulitplixed experiments
##'
##' @param ecol A `character(1)` indicating which column in
##'     `reportData` contains the quantitative information. Default is
##'     `"MS1.Area"`.
##'
##' @param multiplexing A `character(1)` indicating the type of
##'     multiplexing used in the experiment. Provide `"none"` if the
##'     experiment is label-free (default). Available options are:
##'     `"mTRAQ"`.
##'
##' @param ... Further arguments passed to [readQFeatures()].
##'
##' @return An instance of class `QFeatures`. The expression data of
##'     each acquisition run is stored in a separate assay as a
##'     `SummarizedExperiment` object.
##'
##' - The `QFeatures` (see [QFeatures()]) class to read about how to
##'   manipulate the resulting `QFeatures` object.
##'
##' - The [readQFeatures()] function to import other quantitative
##'   data.
##'
##' @author Laurent Gatto, Christophe Vanderaa
##'
##' @importFrom tidyr pivot_wider
##'
##' @export
##'
##' @examples
##'
##' x <- read.delim(MsDataHub::benchmarkingDIA.tsv())
##' ## fix file names
##' x[[1]] <- sub("^.+raw-data\\\\", "", x[[1]])
##' cd <- data.frame(File.Name = unique(x[[1]]))
##' readQFeaturesFromDIANN(colData = cd, reportData = x, ecol = "Ms1.Area")
readQFeaturesFromDIANN <- function(colData, reportData, extractedData = NULL,
                                   ecol = "MS1.Area",
                                   multiplexing = "none", # "none" or "mTRAQ"
                                   ...) {
    diannReportCols <- c("File.Name", "Precursor.Id", "Modified.Sequence")
    if (!all(diannReportCols %in% colnames(reportData)))
        stop("'reportData' is not an expected DIA-NN report table ",
             "output. This function expects the main output file as ",
             "described here: https://github.com/vdemichev/DiaNN#main-output-reference")
    if (!ecol %in% colnames(reportData))
        stop("'", ecol, "' not found in 'reportData'")
    if (!"File.Name" %in% colnames(colData))
        stop("'colData' must contain a column named 'File.Name' that provides ",
             "a link to the 'File.Name' column in 'reportData'")
    if (multiplexing == "none" && !is.null(extractedData))
        stop("Providing 'extractedData' for label-free experiments ",
             "('multiplexed == \"none\"') is not expected. Raise an ",
             "issue if you need this feature: ",
             "https://github.com/UCLouvain-CBIO/scp/issues/new/choose")

    args <- list(...)
    ## Get the label used for the reportData
    if (multiplexing == "mTRAQ") {
        ## Extracted the mTRAQ label from the modified sequence
        reportData$Label <- sub("^.*[Q-](\\d).*$", "\\1", reportData$Modified.Sequence)
        reportData$Precursor.Id <- gsub("\\(mTRAQ.*?\\)", "(mTRAQ)", reportData$Precursor.Id)
        args$sep <- "."
        ## Make sure the colData has the Label column
        if (!"Label" %in% colnames(colData))
            stop("'colData' must contain a column named 'Label' that ",
                 "provides the mTRAQ reagent used to label the ",
                 "samples and/or single cells.")
        if (any(mis <- !colData$Label %in% reportData$Label)) {
            stop("Some labels from 'colData$Label' were not found as",
                 "part of the mTRAQ labels found in ",
                 "'reportData$Modified.Sequence': ",
                 paste0(unique(colData$Label[mis]), collapse = ", "))
        }
        ## Identify which variables are correlated with the run-specific
        ## precursor IDs
        nIds <- length(unique(paste0(reportData$Precursor.Id, reportData$File.Name)))
        nLevels <- sapply(colnames(reportData), function(x) {
            nrow(unique(reportData[, c("Precursor.Id", "File.Name", x)]))
        })
        idCols <- names(nLevels)[nLevels == nIds]
        ## Transform the reportData to a wide format with respect to label
        reportData <- pivot_wider(reportData, id_cols = all_of(idCols),
                                  names_from = "Label",
                                  values_from = ecol)
    } else if (multiplexing == "none") {
        colData$Label <- ecol
        args$sep <- ""
        args$suffix <- ""
    } else {
        stop("The '", multiplexing, "' multiplexing strategy is not ",
             "implemented. Raise an issue if you need this feature: ",
             "https://github.com/UCLouvain-CBIO/scp/issues/new/choose")
    }

    ## Read using readSCP
    out <- do.call(readQFeatures, c(args, list(assayData = reportData,
                                               colAnnotation = colData,
                                               batchCol = "File.Name",
                                               channelCol = "Label")))

    ## Optionally, add the extractedData
    if (!is.null(extractedData)) {
        labs <- unique(colData$Label)
        ## DIA-NN appends the label to the run name
        quantCols <- grep(paste0("[", paste0(labs, collapse = ""), "]$"),
                          colnames(extractedData))
        extractedData <- readSummarizedExperiment(extractedData,
                                                  ecol = quantCols,
                                                  fnames = "Precursor.Id")
        ## Make sure extractedData has the sames samples as reportData
        cnames <- unique(unlist(colnames(out)))
        if (any(mis <- !cnames %in% colnames(extractedData)))
            stop("Some columns present in reportData are not found in ",
                 "extracted data", paste0(cnames[mis], collapse = ", "),
                 "\nAre you sure the two tables were generated from ",
                 "the same experiment?")
        extractedData <- extractedData[, cnames]
        ## Add the assay to the QFeatures object
        anames <- names(out)
        out <- addAssay(out, extractedData, name = "Ms1Extracted")
        out <- addAssayLink(out,
                            from = anames, to = "Ms1Extracted",
                            varFrom = rep("Precursor.Id", length(anames)),
                            varTo = "Precursor.Id")
    }
    out
}
