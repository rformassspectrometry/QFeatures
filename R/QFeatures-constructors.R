##' @title QFeatures from tabular data
##'
##' @description
##'
##' These functions convert tabular data into dedicated data
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
##'   the input table. This multi-set case should be used when the
##'   input table contains data for multiple runs/batches?
##'
##' @param assayData A `data.frame`, or any object that can be coerced
##'     to a `data.frame`, holding the quantitative assay. For
##'     `readSummarizedExperiment()`, this can also be a
##'     `character(1)` to a filename.
##'
##' @param colAnnotation The type of this parameter will define
##'     whether the resulting `QFeatures` object will contain a single
##'     or multiple sets.
##'
##'     For the single-set case, a `numeric` indicating the indices of
##'     the columns to be used as expression values, or a `character`
##'     indicating the names of the columns, or a `logical` indicating
##'     the quantitative assay's columns.
##'
##'     For the multi-set case, a `data.frame` or any object that can
##'     be coerced to a `data.frame`. It is expected to contain all
##'     the sample annotations. Required fields are the acquisition
##'     batch (given by `runCol`) and the acquisition channel within
##'     the batch (e.g. TMT channel, given by
##'     `channelCol`). Additional fields (e.g. sample type,
##'     acquisition date,...) are allowed and will be stored as sample
##'     metadata in the `QFeatures`'s colData slot.
##'
##' @param quantCols A `numeric()`, `logical()` or `character()`
##'     defining the columns of the `assayData` that contain the
##'     quantitative data. Can also be defined in `colAnnotation`.
##'
##' @param runCol For the multi-set case, a `numeric(1)` or
##'     `character(1)` pointing to the column of `assayData` and
##'     `colAnnotation` that contain the batch names. Make sure that
##'     the column name in both tables are either identical and
##'     syntactically valid (if you supply a `character`) or have the
##'     same index (if you supply a `numeric`). Note that characters
##'     are converted to syntactically valid names using `make.names`
##'
##' @param fnames For the single- and multi-set cases, an optional
##'     `character(1)` or `numeric(1)` indicating the column to be
##'     used as feature names.  Note that rownames must be unique in
##'     `QFeatures` sets.
##'
##' @param name For the single-set case, an optional `character(1)` to
##'     name the set in the `QFeatures` object. Default is `quants`.
##'
##' @param removeEmptyCols A `logical(1)`. If true, the function will
##'     remove in each batch the columns that contain only missing
##'     values.
##'
##' @param verbose A `logical(1)` indicating whether the progress of
##'     the data reading and formatting should be printed to the
##'     console. Default is `TRUE`.
##'
##' @return An instance of class `QFeatures` or
##'     [SummarizedExperiment::SummarizedExperiment()]. For the
##'     former, the quantitative sets of each run are stored in
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
##' - The [readQFeaturesFromDIANN()] function to import DIA-NN
##'   quantitative data.
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
##' ## Create a QFeatures object with a single psms set
##' qf1 <- readQFeatures(hlpsms, quantCols = 1:10, name = "psms")
##' qf1
##'
##' ######################################
##' ## Single-set case using a data.frame
##'
##' ## All PSMs were acquired in the same acquisition
##' hlpsms$file <- "File1"
##' hlpsms[1:10, c(1, 2, 10:11, 14, 17, 29)]
##' (colann <- data.frame(file = rep("File1", 10),
##'                       var = rnorm(10),
##'                       quantCols = names(hlpsms)[1:10]))
##' qf2 <- readQFeatures(hlpsms, colAnnotation = colann)
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
##' (colann <- data.frame(runCol = rep(paste0("File", 1:3), each = 10),
##'                       quantCols = names(hlpsms)[1:10]))
##'
##' qf3 <- readQFeatures(hlpsms, colAnnotation = colann, runCol = "file")
##' qf3
##' colData(qf3)
NULL

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

##' @export
##'
##' @rdname readQFeatures
##'
##' @importFrom utils read.csv
##'
##' @param ecol Same as `quantCols` for the single-set case. Available
##'     for backwards compatibility. Default is `NULL`. If both `ecol`
##'     and `colAnnotation` are set, an error is thrown.
##'
##' @param ... Further arguments that can be passed on to [read.csv()]
##'     except `stringsAsFactors`, which is always `FALSE`. Only
##'     applicable to `readSummarizedExperiment()`.
readSummarizedExperiment <- function(assayData,
                                     quantCols = NULL,
                                     fnames = NULL,
                                     ecol = NULL, ...) {
    if (!is.null(ecol)) {
        warning("'ecol' is deprecated, use 'quantCols' instead.")
        if (is.null(quantCols)) quantCols <- ecol
    }
    if (!is.vector(quantCols) || is.list(quantCols))
        stop("'quantCols' must be an atomics vector.")
    if (is.data.frame(assayData)) xx <- assayData
    else {
        args <- list(...)
        args$file <- assayData
        if ("rownames" %in% names(args)) {
            if (is.null(fnames)) fnames <- args$rownames
            args$rownames <- NULL
        }
        args$stringsAsFactors <- FALSE
        xx <- do.call(read.csv, args)
    }
    if (is.character(quantCols) || is.factor(quantCols)) {
        mis <- !quantCols %in% colnames(xx)
        if (any(mis))
            stop("Column identifiers ",
                 paste(quantCols[mis], collapse = ", "),
                 " not recognised among\n",
                 paste(colnames(xx), paste = ", "))
        quantCols <- which(colnames(xx) %in% quantCols)
    } else if (is.logical(quantCols)) {
        if (length(quantCols) != length(xx))
            stop("Length of 'quantCols' and 'assayData' do not match.")
        quantCols <- which(quantCols)
    }
    assay <- as.matrix(xx[, quantCols, drop = FALSE])
    fdata <- DataFrame(xx[, -quantCols, drop = FALSE])

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
##' @rdname readQFeatures
readQFeatures <- function(assayData,
                          colAnnotation = NULL,
                          quantCols = NULL,
                          runCol = NULL,
                          name = "quants",
                          removeEmptyCols = FALSE,
                          verbose = TRUE,
                          ...) {
    if (verbose) message("Checking arguments.")
    assayData <- as.data.frame(assayData)
    if (!is.null(colAnnotation))
        colAnnotation <- data.frame(colAnnotation)
    quantCols <- .checkQuantCols(assayData, colAnnotation, quantCols)
    runs <- .checkRunCol(assayData, colAnnotation, runCol)
    if (verbose) message("Loading data as a 'SummarizedExperiment' object.")
    se <- readSummarizedExperiment(assayData, quantCols, ...)
    rownames(se) <- make.unique(rownames(se))
    if (length(runs)) {
        if (verbose) message("Splitting data in runs.")
        el <- .splitSE(se, runs)
        el <- .createUniqueColnames(el, quantCols)
    } else {
        el <- structure(list(se), .Names = name[1])
    }
    if (removeEmptyCols) el <- .removeEmptyColumns(el)
    if (verbose) message("Formatting sample annotations (colData).")
    colData <- .formatColData(el, colAnnotation, runs, quantCols)
    if (verbose) message("Formatting data as a 'QFeatures' object.")
    QFeatures(experiments = el, colData = colData)
}

## This function will check the quantitation variable inputs. At the
## end, it will return a valid character quantCols, either as provided
## directly by the user, or generated from colAnnotation.
.checkQuantCols <- function(assayData, colAnnotation, quantCols) {
    ## Fail early if both a missing
    if (is.null(colAnnotation) & is.null(quantCols))
        stop("Provide one of 'colAnnotation' or 'quantCols', both mustn't be NULL.")
    ## If we have a colAnnotation data.frame, it must contain a
    ## quantCols column and no quantCols should be provided.
    if (!is.null(colAnnotation)) {
        if (!"quantCols" %in% colnames(colAnnotation) &&
            length(quantCols) > 1)
            stop("'colAnnotation' must contain a column called 'quantCols'")
    }
    if (is.null(quantCols)) {
        ## if (is.null(colAnnotation))
        ##     stop("'quantCols' and 'colAnnotation' cannot both be NULL.")
        if (!"quantCols" %in% colnames(colAnnotation))
            stop("When 'quantCols' is NULL, 'colAnnotation' must ",
                 "contain a column called 'quantCols'.")
        quantCols <- unique(colAnnotation$quantCols)
    }
    if (is.numeric(quantCols) || is.logical(quantCols))
        quantCols <- colnames(assayData)[quantCols]
    mis <- quantCols[!quantCols %in% colnames(assayData)]
    if (length(mis))
        stop("Some column names in 'quantCols' are not found ",
             "in 'assayData': ", paste0(mis, collapse = ", "), ".")
    quantCols
}


## This function will check the batch/run variable inputs. At the end,
## it will return a vector of runs/batches or NULL.
.checkRunCol <- function(assayData, colAnnotation, runCol) {
    if (is.null(runCol)) return(NULL)
    if (length(runCol) > 1)
        stop("'runCol' is a vector. Please provide the name of a ",
             "column in 'assayData'.")
    if (!runCol %in% colnames(assayData))
        stop("'", runCol, "' (provided as 'runCol') not found ",
             "in 'assayData'.")
    runs <- assayData[[runCol]]
    if (!is.null(colAnnotation)) {
        if (!"runCol" %in% colnames(colAnnotation))
            stop("When 'runCol' is not NULL, 'colAnnotation' must ",
                 "contain a column called 'runCol'.")
        mis <- !runs %in% colAnnotation$runCol
        if (any(mis)) {
            warning("Some runs are missing in 'colAnnot': ",
                    paste0(unique(runs[mis]), collapse = ", "))
        }
    }
    assayData[[runCol]]
}

##' Split SummarizedExperiment into an ExperimentList
##'
##' The fonction creates an [ExperimentList] containing
##' [SummarizedExperiment] objects from a [SummarizedExperiment]
##' object (also works with [SingleCellExperiment] objects). `f` is
##' used to split `x`` along the rows (`f`` was a feature variable
##' name) or samples/columns (f was a phenotypic variable name). If f
##' is passed as a factor, its length will be matched to nrow(x) or
##' ncol(x) (in that order) to determine if x will be split along the
##' features (rows) or sample (columns). Hence, the length of f must
##' match exactly to either dimension.
##'
##' This function is not exported and was initially available as
##' scp::.splitSCE().
##'
##' @param x a single [SummarizedExperiment] object
##'
##' @param f a factor or a character of length 1. In the latter case,
##'     `f` will be matched to the row and column data variable names
##'     (in that order). If a match is found, the respective variable
##'     is extracted, converted to a factor if needed.
##' @noRd
.splitSE <- function(x, f) {
    ## Check that f is a factor
    if (length(f) == 1) {
        if (f %in% colnames(rowData(x))) {
            f <- rowData(x)[, f]
        }
        else if (f %in% colnames(colData(x))) {
            f <- colData(x)[, f]
        }
        else {
            stop("'", f, "' not found in rowData or colData")
        }
    }
    ## Check that the factor matches one of the dimensions
    if (!length(f) %in% dim(x))
        stop("length(f) not compatible with dim(x).")
    if (length(f) == nrow(x)) { ## Split along rows
        xl <- lapply(split(rownames(x), f = f), function(i) x[i, ])
    } else { ## Split along columns
        xl <- lapply(split(colnames(x), f = f), function(i) x[, i])
    }
    ## Convert list to an ExperimentList
    do.call(ExperimentList, xl)
}

.createUniqueColnames <- function(el, quantCols) {
    if (length(quantCols) == 1)  suffix <- ""
    else suffix <- paste0("_", quantCols)
    for (i in seq_along(el)) {
        colnames(el[[i]]) <- paste0(names(el)[[i]], suffix)
    }
    el
}

.removeEmptyColumns <- function(el) {
    for (i in seq_along(el)) {
        sel <- colSums(is.na(assay(el[[i]]))) != nrow(el[[i]])
        el[[i]] <- el[[i]][, sel]
    }
    el
}

## This function will create a colData from the different (possibly
## missing, i.e. NULL) arguments
.formatColData <- function(el, colAnnotation, runs, quantCols) {
    sampleNames <- unlist(lapply(el, colnames), use.names = FALSE)
    if (is.null(colAnnotation))
        return(DataFrame(row.names = sampleNames))
    if (!length(runs)) {
        ## rownames(colAnnotation) <- colAnnotation$quantCols
        rownames(colAnnotation) <- sampleNames
    } else {
        if (length(quantCols) == 1) {
            rownames(colAnnotation) <- colAnnotation$runCol
        } else {
            rownames(colAnnotation) <- paste0(colAnnotation$runCol, "_", colAnnotation$quantCols)
        }
    }
    colData <- colAnnotation[sampleNames, , drop = FALSE]
    rownames(colData) <- sampleNames ## clean NA in rownames
    colData
}

##' @title  Read DIA-NN output as a QFeatures objects
##'
##' @description
##'
##' This function takes the output tables from DIA-NN and converts them
##' into a `QFeatures` object.
##'
##' @param assayData A `data.frame` or any object that can be coerced
##'     to a data.frame that contains the data from the `Report.tsv`
##'     file generated by DIA-NN.
##'
##' @param colAnnotation A `data.frame` or any object that can be
##'     coerced to a `data.frame`. `colAnnotation` is expected to
##'     contains all the sample annotations. We require the table to
##'     contain a column called `File.Name` that links to the
##'     `File.Name` in the DIA-NN report table. If `multiplexing =
##'     "mTRAQ"`, we require a second column called `Label` that links
##'     the label to the sample (the labels identified by DIA-NN can
##'     be retrieved from `Modified Sequence` column in the report
##'     table).
##'
##' @param extractedData A data.frame or any object that can be coerced
##'     to a data.frame that contains the data from the `*_ms1_extracted.tsv`
##'     file generated by DIA-NN. This argument is optional and is
##'     currently only applicable for mTRAQ multiplexed experiments
##'     where DIA-NN was run using the `plexdia` module.
##'
##' @param quantCols A `character(1)` indicating which column in
##'     `assayData` contains the quantitative information. Default is
##'     `"Ms1.Area"`.
##'
##' @param multiplexing A `character(1)` indicating the type of
##'     multiplexing used in the experiment. Provide `"none"` if the
##'     experiment is label-free (default). Alternative options are:
##'     `"mTRAQ"`.
##'
##' @param ecol Same as `quantCols` for the single-set case. Available
##'     for backwards compatibility. Default is `NULL`. If both `ecol`
##'     and `colAnnotation` are set, an error is thrown.
##'
##' @param ... Further arguments passed to [readQFeatures()].
##'
##' @return An instance of class `QFeatures`. The quantiative data of
##'     each acquisition run is stored in a separate set as a
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
##' @importFrom tidyselect all_of
##'
##' @export
##'
##' @examples
##'
##' ## x <- read.delim(MsDataHub::benchmarkingDIA.tsv())
##' ## fix file names
##' ## x[[1]] <- sub("^.+raw-data\\\\", "", x[[1]])
##' ## cd <- data.frame(File.Name = unique(x[[1]]))
##' ## readQFeaturesFromDIANN(colAnnotation = cd, assayData = x)
readQFeaturesFromDIANN <- function(assayData, colAnnotation, extractedData = NULL,
                                   quantCols = "Ms1.Area", multiplexing = "none",
                                   ecol = NULL,
                                   ...) {
    suppArgs <- .checkDiannArguments(
        colAnnotation, assayData, extractedData, ecol, multiplexing, ...
    )
    if (multiplexing == "mTRAQ") {
        assayData <- .formatMtraqReportData(assayData, colAnnotation, ecol)
    } else if (multiplexing == "none") {
        colAnnotation$Label <- ecol
    }
    allArgs <- c(suppArgs, list(
        assayData = assayData, colAnnotation = colAnnotation,
        runCol = "File.Name", channelCol = "Label"
    ))
    out <- do.call(readQFeatures, allArgs)
    if (!is.null(extractedData)) {
        out <- .addDiannExtractedData(out, extractedData)
    }
    out
}

## Internal function that checks whether the provided arguments match
## the expected input that is generated by DIA-NN.
## Parameter description is the same as for `readSCPfromDIANN()`
.checkDiannArguments <- function(colAnnotation, assayData, extractedData,
                                 ecol, multiplexing, ...) {
    diannReportCols <- c("File.Name", "Precursor.Id", "Modified.Sequence")
    if (!all(diannReportCols %in% colnames(assayData)))
        stop("'assayData' is not an expected DIA-NN report table ",
             "output. This function expects the main output file as ",
             "described here: https://github.com/vdemichev/DiaNN#main-output-reference")
    if (!ecol %in% colnames(assayData))
        stop("'", ecol, "' not found in 'assayData'")
    if (!"File.Name" %in% colnames(colAnnotation))
        stop("'colAnnotation' must contain a column named 'File.Name' that provides ",
             "a link to the 'File.Name' column in 'assayData'")
    if (multiplexing == "none" && !is.null(extractedData))
        stop("Providing 'extractedData' for label-free experiments ",
             "('multiplexed == \"none\"') is not expected. Raise an ",
             "issue if you need this feature: ",
             "https://github.com/UCLouvain-CBIO/scp/issues/new/choose")
    .checkDiannArgumentsDots(multiplexing, ...)
}

## Internal function that adapts the dots arguments (that will be used
## by `readSCP()`) depending on the multiplexing approach used.
.checkDiannArgumentsDots <- function(multiplexing, ...) {
    suppArgs <- list(...)
    if (multiplexing == "mTRAQ") {
        suppArgs$sep <- "."
    } else if (multiplexing == "none") {
        suppArgs$sep <- ""
        suppArgs$suffix <- ""
    } else {
        stop("The '", multiplexing, "' multiplexing strategy is not ",
             "implemented. Raise an issue if you need this feature: ",
             "https://github.com/UCLouvain-CBIO/scp/issues/new/choose")
    }
    suppArgs
}

## (Only for mTRAQ multiplexing!) Internal function that extracts the
## mTRAQlabels from the peptide sequence, removes the mTRAQ annotation
## from the precursor ID, identifies constant columns within precursor
## and puts the quantification data for different mTRAQ labels in
## separate columns (wide format).
.formatMtraqReportData <- function(assayData, colAnnotation, ecol) {
    assayData$Label <-
        sub("^.*[Q-](\\d).*$", "\\1", assayData$Modified.Sequence)
    assayData$Precursor.Id <-
        gsub("\\(mTRAQ.*?\\)", "(mTRAQ)", assayData$Precursor.Id)
    .checkLabelsInColData(colAnnotation, assayData)
    idCols <- .findPrecursorVariables(assayData)
    pivot_wider(
        assayData, id_cols = all_of(idCols),
        names_from = "Label", values_from = all_of(ecol)
    )
}

## Internal function that identifies which variables in the report
## data are constant within each precursor (with each run).
.findPrecursorVariables <- function(assayData) {
    precIds <- paste0(assayData$Precursor.Id, assayData$File.Name)
    nUniqueIds <- length(unique(precIds))
    nLevels <- sapply(colnames(assayData), function(x) {
        nrow(unique(assayData[, c("Precursor.Id", "File.Name", x)]))
    })
    names(nLevels)[nLevels == nUniqueIds]
}

## Internal function that ensures that the assayData and the
## colAnnotation are correctly linked.
.checkLabelsInColData <- function(colAnnotation, assayData) {
    if (!"Label" %in% colnames(colAnnotation))
        stop("'colAnnotation' must contain a column named 'Label' that ",
             "provides the mTRAQ reagent used to label the ",
             "samples and/or single cells.")
    if (any(mis <- !colAnnotation$Label %in% assayData$Label)) {
        stop("Some labels from 'colAnnotation$Label' were not found as",
             "part of the mTRAQ labels found in ",
             "'assayData$Modified.Sequence': ",
             paste0(unique(colAnnotation$Label[mis]), collapse = ", "))
    }
    NULL
}

## Internal function that adds the extractedData to a QFeatures
## object. The functions first converts the extractedData to a
## SummarizedExperiment objects and subsets the SCE for the set of
## shared samples. The added assay is automatically linked (using
## AssayLinks) to the assayData.
## Developer's note: the function assumes that DIA-NN creates sample
## names in the extracted data by appending the labels to the run
## names
## @param object A QFeatures object containing DIA-NN report data, as
##     generated by readQFeatures.
## @param extractedData A data.frame or any object that can be coerced
##     to a data.frame that contains the data from the
##     `*_ms1_extracted.tsv` file generated by DIA-NN.
.addDiannExtractedData <- function(object, extractedData) {
    quantColPattern <- paste0(unique(object$Label), "$", collapse = "|")
    quantCols <- grep(quantColPattern, colnames(extractedData))
    extractedData <- readSummarizedExperiment(
        extractedData, ecol = quantCols, fnames = "Precursor.Id"
    )
    extractedData <- .keepSharedSamples(extractedData, object)
    object <- addAssay(object, extractedData, name = "Ms1Extracted")
    addAssayLink(
        object,
        from = grep("Ms1Extracted", names(object), invert = TRUE),
        to = "Ms1Extracted",
        varFrom = rep("Precursor.Id", length(names(object)) - 1),
        varTo = "Precursor.Id"
    )
}

## Internal functions that subsets the extractedData to keep only
.keepSharedSamples <- function(extractedData, object) {
    cnames <- unique(unlist(colnames(object)))
    if (any(mis <- !cnames %in% colnames(extractedData)))
        stop("Some columns present in assayData are not found in ",
             "extracted data", paste0(cnames[mis], collapse = ", "),
             "\nAre you sure the two tables were generated from ",
             "the same experiment?")
    extractedData[, cnames]
}



## @param channelCol For the multi-set case, a `numeric(1)` or
##     `character(1)` pointing to the column of `colAnnotation` that
##     contains the column names of the quantitative data in
##     `assayData` (see example). Useful for multiplexed experiments
##     such as mTRAQ or TMT.
