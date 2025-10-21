##' @title QFeatures from tabular data
##'
##' @description
##'
##' These functions convert tabular data into dedicated data
##' objets. The [readSummarizedExperiment()] function takes a file
##' name or `data.frame` and converts it into a
##' [SummarizedExperiment()] object.  The [readQFeatures()] function
##' takes a `data.frame` and converts it into a `QFeatures` object
##' (see [QFeatures()] for details). For the latter, two use-cases
##' exist:
##'
##' - The single-set case will generate a `QFeatures` object with a
##'   single `SummarizedExperiment` containing all features of the
##'   input table.
##'
##' - The multi-set case will generate a `QFeatures` object containing
##'   multiple `SummarizedExperiment`s, resulting from splitting the
##'   input table. This multi-set case is generally used when the
##'   input table contains data from multiple runs/batches.
##'
##' @details
##'
##' The single- and multi-set cases are defined by the `quantCols` and
##' `runCol` parameters, whether passed by the `quantCols` and
##' `runCol` vectors and/or the `colData` `data.frame` (see below).
##'
##' ## Single-set case
##'
##' The quantitative data variables are defined by the `quantCols`.
##' The single-set case can be represented schematically as shown
##' below.
##'
##' ```
##' |------+----------------+-----------|
##' | cols | quantCols 1..N | more cols |
##' | .    | ...            | ...       |
##' | .    | ...            | ...       |
##' | .    | ...            | ...       |
##' |------+----------------+-----------|
##' ```
##'
##' Note that every `quantCols` column contains data for a single
##' sample. The single-set case is defined by the absence of any
##' `runCol` input (see next section). We here provide a
##' (non-exhaustive) list of typical data sets that fall under the
##' single-set case:
##'
##' - Peptide- or protein-level label-free data (bulk or single-cell).
##' - Peptide- or protein-level multiplexed (e.g. TMT) data (bulk or
##'   single-cell).
##' - PSM-level multiplexed data acquired in a single MS run (bulk or
##'   single-cell).
##' - PSM-level data from fractionation experiments, where each
##'   fraction of the same sample was acquired with the same
##'   multiplexing label.
##'
##' ## Multi-set case
##'
##' A run/batch variable, `runCol`, is required to import multi-set
##' data. The multi-set case can be represented schematically as shown
##' below.
##'
##' ```
##' |--------+------+----------------+-----------|
##' | runCol | cols | quantCols 1..N | more cols |
##' |   1    | .    | ...            | ...       |
##' |   1    | .    | ...            | ...       |
##' |--------+------+----------------+-----------|
##' |   2    | .    | ...            | ...       |
##' |--------+------+----------------+-----------|
##' |   .    | .    | ...            | ...       |
##' |--------+------+----------------+-----------|
##' ```
##'
##' Every `quantCols` column contains data for multiple samples
##' acquired in different runs. The multi-set case applies when
##' `runCol` is provided, which will determine how the table is split
##' into multiple sets.
##'
##' We here provide a (non-exhaustive) list of typical data sets that
##' fall under the multi-set case:
##'
##' - PSM- or precursor-level multiplexed data acquired in multiple
##'   runs (bulk or single-cell)
##' - PSM- or precursor-level label-free data acquired in multiple
##'   runs (bulk or single-cell)
##' - DIA-NN data (see also [readQFeaturesFromDIANN()]).
##'
##' ## Adding sample annotations with `colData`
##'
##' We recommend providing sample annotations when creating a
##' `QFeatures` object. The `colData` is a table in which each row
##' corresponds to a sample and each column provides information about
##' the samples. There is no restriction on the number of columns and
##' on the type of data they should contain. However, we impose one or
##' two columns (depending on the use case) that allow to link the
##' annotations of each sample to its quantitative data:
##'
##' - Single-set case: the `colData` must contain a column named
##'   `quantCols` that provides the names of the columns in
##'   `assayData` containing quantitative values for each sample (see
##'   single-set cases in the examples).
##'
##' - Multi-set case: the `colData` must contain a column named
##'   `quantCols` that provides the names of the columns in
##'   `assayData` with the quantitative values for each sample, and a
##'   column named `runCol` that provides the MS runs/batches in which
##'   each sample has been acquired. The entries in
##'   `colData[["runCol"]]` are matched against the entries provided
##'   by `assayData[[runCol]]`.
##'
##' When the `quantCols` argument is not provided to
##' `readQFeatures()`, the function will automatically determine the
##' `quantCols` from `colData[["quantCols"]]`. Therefore, `quantCols`
##' and `colData` cannot be both missing.
##'
##' Samples that are present in `assayData` but absent
##' `colData` will lead to a warning, and the missing entries will be
##' automatically added to the `colData` and filled with `NA`s.
##'
##' When using the `quantCols` and `runCol` arguments only
##' (without `colData`), the `colData` contains zero
##' columns/variables.
##'
##' @param assayData A `data.frame`, or any object that can be coerced
##'     into a `data.frame`, holding the quantitative assay. For
##'     `readSummarizedExperiment()`, this can also be a
##'     `character(1)` pointing to a filename. This `data.frame` is
##'     typically generated by an identification and quantification
##'     software, such as Sage, Proteome Discoverer, MaxQuant, ...
##'
##' @param colData A `data.frame` (or any object that can be coerced
##'     to a `data.frame`) containing sample/column annotations,
##'     including `quantCols` and `runCol` (see details).
##'
##' @param quantCols A `numeric()`, `logical()` or `character()`
##'     defining the columns of the `assayData` that contain the
##'     quantitative data. This information can also be defined in
##'     `colData` (see details).
##'
##' @param runCol For the multi-set case, a `numeric(1)` or
##'     `character(1)` pointing to the column of `assayData` (and
##'     `colData`, is set) that contains the runs/batches. Make sure
##'     that the column name in both tables are identical and
##'     syntactically valid (if you supply a `character`) or have the
##'     same index (if you supply a `numeric`). Note that characters
##'     are converted to syntactically valid names using `make.names`
##'
##' @param fnames For the single- and multi-set cases, an optional
##'     `character(1)` or `numeric(1)` indicating the column to be
##'     used as feature names.  Note that rownames must be unique
##'     within `QFeatures` sets. Default is `NULL`.
##'
##' @param name For the single-set case, an optional `character(1)` to
##'     name the set in the `QFeatures` object. Default is `quants`.
##'
##' @param removeEmptyCols A `logical(1)`. If `TRUE`, quantitative
##'     columns that contain only missing values are removed.
##'
##' @param verbose A `logical(1)` indicating whether the progress of
##'     the data reading and formatting should be printed to the
##'     console. Default is `TRUE`.
##'
##' @param ecol Same as `quantCols`. Available for backwards
##'     compatibility. Default is `NULL`. If both `ecol` and `colData`
##'     are set, an error is thrown.
##'
##' @param ... Additional parameters passed to
##'     `readSummarizedExperiment()` by `readQFeatures()` and
##'     [read.csv()] by `readSummarizedExperiment()`.
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
##' colData(qf1)
##'
##' ######################################
##' ## Single-set case with colData.
##'
##' (coldat <- data.frame(var = rnorm(10),
##'                       quantCols = names(hlpsms)[1:10]))
##' qf2 <- readQFeatures(hlpsms, colData = coldat)
##' qf2
##' colData(qf2)
##'
##' ######################################
##' ## Multi-set case.
##'
##' ## Let's simulate 3 different files/batches for that same input
##' ## data.frame, and define a colData data.frame.
##'
##' hlpsms$file <- paste0("File", sample(1:3, nrow(hlpsms), replace = TRUE))
##' hlpsms[1:10, c(1, 2, 10:11, 14, 17, 29)]
##'
##' qf3 <- readQFeatures(hlpsms, quantCols = 1:10, runCol = "file")
##' qf3
##' colData(qf3)
##'
##'
##' ######################################
##' ## Multi-set case with colData.
##'
##' (coldat <- data.frame(runCol = rep(paste0("File", 1:3), each = 10),
##'                       var = rnorm(10),
##'                       quantCols = names(hlpsms)[1:10]))
##' qf4 <- readQFeatures(hlpsms, colData = coldat, runCol = "file")
##' qf4
##' colData(qf4)
NULL


##' @export
##'
##' @rdname readQFeatures
##'
##' @importFrom utils read.csv
##'
##' @param ... Further arguments that can be passed on to [read.csv()]
##'     except `stringsAsFactors`, which is always `FALSE`. Only
##'     applicable to `readSummarizedExperiment()`.
readSummarizedExperiment <- function(assayData,
                                     quantCols = NULL,
                                     fnames = NULL,
                                     ecol = NULL, ...) {
    quantCols <- .checkWarnEcol(quantCols, ecol)
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
                          colData = NULL,
                          quantCols = NULL,
                          runCol = NULL,
                          name = "quants",
                          removeEmptyCols = FALSE,
                          verbose = TRUE,
                          ecol = NULL,
                          fnames = NULL,
                          ...) {
    if (verbose) message("Checking arguments.")
    assayData <- as.data.frame(assayData)
    if (!is.null(colData))
        colData <- data.frame(colData)
    quantCols <- .checkWarnEcol(quantCols, ecol)
    quantCols <- .checkQuantCols(assayData, colData, quantCols)
    runs <- .checkRunCol(assayData, colData, runCol)
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
    colData <- .formatColData(el, colData, runs)
    if (verbose) message("Formatting data as a 'QFeatures' object.")
    ans <- QFeatures(experiments = el, colData = colData)
    if (!is.null(fnames)) {
        if (verbose) message("Setting assay rownames.")
        ans <- .setAssayRownames(ans, fnames)
    }
    setQFeaturesType(ans, "bulk")
}


## ecol will be deprecated next release. This function warns if ecol
## is used (i.e. is not NULL), then sets quantCols with the value of
## ecol if quantCols wasn't used.
.checkWarnEcol <- function(quantCols, ecol) {
    if (!is.null(ecol)) {
        if (!is.null(quantCols))
            stop("'quantCols' and 'ecols' can't be defined together. ",
                 "Use 'quantCols' only.")
        warning("'ecol' is deprecated, use 'quantCols' instead.")
        if (is.null(quantCols))
            quantCols <- ecol
    }
    quantCols
}

## This function will check the quantitation variable inputs. At the
## end, it will return a valid character quantCols, either as provided
## directly by the user, or generated from colData.
.checkQuantCols <- function(assayData, colData, quantCols) {
    ## Fail early if both a missing
    if (is.null(colData) & is.null(quantCols))
        stop("Provide one of 'colData' or 'quantCols', both mustn't be NULL.")
    ## If we have a colData data.frame, it must contain a
    ## quantCols column and no quantCols should be provided.
    if (!is.null(colData)) {
        if (!"quantCols" %in% colnames(colData) &&
            length(quantCols) > 1)
            stop("'colData' must contain a column called 'quantCols'")
    }
    if (is.null(quantCols)) {
        ## if (is.null(colData))
        ##     stop("'quantCols' and 'colData' cannot both be NULL.")
        if (!"quantCols" %in% colnames(colData))
            stop("When 'quantCols' is NULL, 'colData' must ",
                 "contain a column called 'quantCols'.")
        quantCols <- unique(colData$quantCols)
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
## it will return a vector of runs/batches or NULL (single-set
## case). Possible inputs combinations are:
##
## - `runCol` is NULL: single-set case
## - `runCol` only (i.e. `colAnnotion` is NULL), of length 1, refering
##   to a variable in `assayData`.
## - `runCol` and `colData`: in this case,
##   `colData$runCol` must exist.
.checkRunCol <- function(assayData, colData, runCol) {
    ## No runCol provided: single-set case
    if (is.null(runCol)) return(NULL)
    ## We have a runCol argument: multi-set case
    if (length(runCol) > 1)
        stop("'runCol' must contain the name of a single column ",
             "in 'assayData'.")
    if (!runCol %in% colnames(assayData))
        stop("'", runCol, "' (provided as 'runCol') not found ",
             "in 'assayData'.")
    runs <- assayData[[runCol]]
    if (!is.null(colData)) {
        ## We have a colData argument
        if (!"runCol" %in% colnames(colData))
            stop("When 'runCol' is not NULL, 'colData' must ",
                 "contain a column called 'runCol'.")
        mis <- !runs %in% colData$runCol
        if (any(mis)) {
            warning("Some runs are missing in 'colData': ",
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
.formatColData <- function(el, colData, runs) {
    sampleNames <- unlist(lapply(el, colnames), use.names = FALSE)
    if (is.null(colData))
        return(DataFrame(row.names = sampleNames))
    # assign colData rownames to match colData to sampleNames
    # colData$quantCols presence was checked by .checkQuantCols
    if (is.null(runs)) {
        # use colData$quantCols as rownames
        rownames(colData) <- colData$quantCols
    } else {
        # run information is present, colData should match individual runs
        # the presence of colData$runCol was checked by .checkRunCol
        if (any(duplicated(colData$runCol))) {
            # quantCols as postfix if runCol is duplicated
            newRownames <- paste0(colData$runCol, "_", colData$quantCols)
            if (any(duplicated(newRownames)))
                stop("There are duplicated samples (runCol-quantCols combinations) in the colData table.")
            rownames(colData) <- newRownames
        } else {
            rownames(colData) <- colData$runCol
        }
    }
    # match colData to sampleNames
    colData <- colData[sampleNames, , drop = FALSE]
    rownames(colData) <- sampleNames ## clean NA in rownames
    colData
}


## This function sets the assay rownames. We use it in readQFeatures() when
## fnames is used. Note that we don't want to export it to avoid messing
## with rownames and assayLinks.
.setAssayRownames <- function(object, fcol) {
    stopifnot(inherits(object, "MultiAssayExperiment"))
    ok <- lapply(rowData(object),
                 function(x) stopifnot(fcol %in% names(x)))
    expl <- lapply(experiments(object),
                   function(x) {
                       rn <- rowData(x)[[fcol]]
                       if (anyDuplicated(rn))
                           rn <- make.unique(rn)
                       rownames(x) <- rn
                       x
                   })
    experiments(object) <- List(expl)
    object
}
