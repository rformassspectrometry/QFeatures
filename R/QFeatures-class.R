##' @title Quantitative MS QFeatures
##'
##' @description
##'
##' Conceptually, a `QFeatures` object holds a set of *assays*, each
##' composed of a `matrix` (or `array`) containing quantitative data
##' and row annotations (meta-data).  The number and the names of the
##' columns (samples) must always be the same across the assays, but
##' the number and the names of the rows (features) can vary. The
##' assays are typically defined as `SummarizedExperiment` objects. In
##' addition, a `QFeatures` object also uses a single `DataFrame` to
##' annotate the samples (columns) represented in all the matrices.
##'
##' The `QFeatures` class extends the
##' [MultiAssayExperiment::MultiAssayExperiment] and inherits all
##' the functionality of the
##' [MultiAssayExperiment::MultiAssayExperiment] class.
##'
##' A typical use case for such `QFeatures` object is to represent
##' quantitative proteomics (or metabolomics) data, where different
##' assays represent quantitation data at the PSM (the main assay),
##' peptide and protein level, and where peptide values are computed
##' from the PSM data, and the protein-level data is calculated based
##' on the peptide-level values. The largest assay (the one with the
##' highest number of features, PSMs in the example above) is
##' considered the main assay.
##'
##' The recommended way to create `QFeatures` objects is the use the
##' `readQFeatures()` function, that creates an instance from tabular
##' data. The `QFeatures` constructor can be used to create objects
##' from their bare parts.  It is the user's responsability to make
##' sure that these match the class validity requirements.
##'
##' @section Constructors:
##'
##' - `QFeatures(..., assayLinks)` allows the manual construction of
##'   objects. It is the user's responsability to make sure these
##'   comply. The arguments in `...` are those documented in
##'   [MultiAssayExperiment::MultiAssayExperiment()]. For details
##'   about `assayLinks`, see [AssayLinks]. An example is shown below.
##'
##' - The [readQFeatures()] function constructs a `QFeatures` object
##'   from text-based spreadsheet or a `data.frame` used to generate
##'   an assay. See the function manual page for details and an
##'   example.
##'
##' @section Accessors:
##'
##' - The `QFeatures` class extends the
##'   [MultiAssayExperiment::MultiAssayExperiment] class and inherits
##'   all its accessors and replacement methods.
##'
##' - The `rowDataNames` accessor returns a list with the `rowData`
##'   variable names.
##'   
##' - The `longFormat` accessor takes a `QFeatures` object and returns 
##'   it in a long format `DataFrame`. Each quantitative value is 
##'   reported on a separate line. `colData` and `rowData` data can 
##'   also be added. This function is an extension of the `longFormat` 
##'   function in the [MultiAssayExperiment::MultiAssayExperiment].
##'
##' - The `getRowData` function extracts the `rowData` from one or more
##'   assays in a `QFeatures` object. The output is a single `DFrame` 
##'   table that contains the `rowData` for all selected assays. On the
##'   other hand, the `setRowData` function modifies columns in the 
##'   `rowData` from one or more assays in `QFeatures` object. This 
##'   function can also be used to remove columns in the `rowData` 
##'   (see the vignette for some examples).
##' 
##' @section Adding assays:
##'
##' - The [aggregateFeatures()] function creates a new assay by
##'   aggregating features of an existing assay.
##'
##' - `addAssay(x, y, name, assayLinks)`: Adds a new assay (or
##'   list of assays) `y` to the `QFeatures` instance `x`. `name`
##'   is a `character(1)` naming the single assay (default is
##'   `"newAssay"`), and is ignored if `y` is a list of
##'   assays. `assayLinks` is an optional [AssayLinks].
##'
##' @section Subsetting:
##'
##' - QFeatures object can be subset using the `x[i, j, k, drop =
##'   TRUE]` paradigm. See the argument descriptions for details.
##'
##' - The [subsetByFeature()] function can be used to subset a
##'   `QFeatures` object using one or multiple feature names that will
##'   be matched across different assays, taking the aggregation
##'   relation between assays.
##'
##' - The `selectRowData(x, rowvars)` function can be used to
##'   select a limited number of `rowData` columns of interest named
##'   in `rowvars` in the `x` instance of class `QFeatures`. 
##' 
##' @section Manipulating assays:
##' 
##' - The `replaceRowDataCols` function replaces one or more columns 
##'   in the `rowData` of one or more assays. Note that if the column
##'   to replace does not exist, a new column will get added to the
##'   `rowData`. 
##' - The `removeRowDataCols` function removes one or more columns 
##'   in the `rowData` of one or more assays. If the column to remove
##'   does not exist in one of the assays, the `rowData` of that assay
##'   will stay unchanged. 
##' 
##' @param i `character()`, `integer()`, `logical()` or `GRanges()`
##'     object for subsetting by rows.
##'
##' @param j `character()`, `logical()`, or `numeric()` vector for
##'     subsetting by `colData` rows.
##'
##' @param k `character()`, `logical()`, or `numeric()` vector for
##'     subsetting by assays
##'
##' @param drop logical (default `TRUE`) whether to drop empty assay
##'     elements in the `ExperimentList`.
##'
##'
##' @return See individual method description for the return value.
##'
##' @seealso
##'
##' - The [readQFeatures()] constructor and the [aggregateFeatures()]
##'   function. The *QFeatures* vignette provides an extended example.
##'
##' - The [QFeatures-filtering] manual page demonstrates how to filter
##'   features based on their rowData.
##'
##' - The [missing-data] manual page to manage missing values in
##'   `QFeatures` objects.
##'
##' - The [QFeatures-processing] and [aggregateFeatures()] manual pages
##'   and *Processing* vignette describe common quantitative data
##'   processing methods using in quantitative proteomics.
##'
##' @import MultiAssayExperiment ProtGenerics
##'
##' @name QFeatures
##'
##' @aliases QFeatures QFeatures-class class:QFeatures addAssay dims,QFeatures-method show,QFeatures-method [,QFeatures,ANY,ANY,ANY-method [,QFeatures,character,ANY,ANY-method
##'
##' @aliases rowDataNames selectRowData
##'
##' @rdname QFeatures-class
##'
##' @exportClass QFeatures
##'
##' @author Laurent Gatto
##'
##' @examples
##' ## ------------------------
##' ## An empty QFeatures object
##' ## ------------------------
##'
##' QFeatures()
##'
##' ## -----------------------------------
##' ## Creating a QFeatures object manually
##' ## -----------------------------------
##'
##' ## two assays (matrices) with matching column names
##' m1 <- matrix(1:40, ncol = 4)
##' m2 <- matrix(1:16, ncol = 4)
##' sample_names <- paste0("S", 1:4)
##' colnames(m1) <- colnames(m2) <- sample_names
##' rownames(m1) <- letters[1:10]
##' rownames(m2) <- letters[1:4]
##'
##' ## two corresponding feature metadata with appropriate row names
##' df1 <- DataFrame(Fa = 1:10, Fb = letters[1:10],
##'                  row.names = rownames(m1))
##' df2 <- DataFrame(row.names = rownames(m2))
##'
##' (se1 <- SummarizedExperiment(m1, df1))
##' (se2 <- SummarizedExperiment(m2, df2))
##'
##' ## Sample annotation (colData)
##' cd <- DataFrame(Var1 = rnorm(4),
##'                 Var2 = LETTERS[1:4],
##'                 row.names = sample_names)
##'
##' el <- list(assay1 = se1, assay2 = se2)
##' fts1 <- QFeatures(el, colData = cd)
##' fts1
##' fts1[[1]]
##' fts1[["assay1"]]
##'
##' ## Rename assay
##' names(fts1) <- c("se1", "se2")
##'
##' ## Add an assay
##' fts1 <- addAssay(fts1, se1[1:2, ], name = "se3")
##' 
##' ## Get the assays feature metadata
##' rowData(fts1) 
##'
##' ## Keep only the Fa variable
##' selectRowData(fts1, rowvars = "Fa")
##' 
##' ## -----------------------------------
##' ## See ?readQFeatures to create a
##' ## QFeatures object from a data.frame
##' ## or spreadsheet.
##' ## -----------------------------------
##' 
NULL


## ----------------------------------
## QFeatures Class ChangeLog
##
## Version 0.1:
##  - Contains MatchedAssayExperiment
## Version 0.2:
##  - Contains MultiAssayExperiment (see issue 46)
## Version 0.3:
##  - Rename to QFeatures (see issue 89)

setClass("QFeatures",
         contains = "MultiAssayExperiment",
         slots = c(version = "character",
                   assayLinks = "AssayLinks"),
         prototype = prototype(
             version = "0.3"))


##' @rdname QFeatures-class
##' @param object An instance of class `QFeatures`.
##' @exportMethod show
setMethod("show", "QFeatures",
          function(object) {
              if (isEmpty(object)) {
                  cat(sprintf("A empty instance of class %s", class(object)), "\n")
                  return(NULL)
              }
              n <- length(object)
              cat(sprintf("An instance of class %s", class(object)), "containing", n, "assays:")
              el <- experiments(object)
              o_class <- class(el)
              elem_cl <- vapply(el, class, character(1L))
              o_len <- length(el)
              o_names <- names(el)
              featdim <- vapply(el, FUN = function(obj) {
                  dim(obj)[1]
              }, FUN.VALUE = integer(1L))
              sampdim <- vapply(el, FUN = function(obj) {
                  dim(obj)[2]
              }, FUN.VALUE = integer(1L))
              if (n <= 7) {
                  cat(sprintf("\n [%i] %s: %s with %s rows and %s columns",
                              seq(o_len), o_names, elem_cl,
                              featdim, sampdim), "\n")
              } else {
                  cat(sprintf("\n [%i] %s: %s with %s rows and %s columns",
                              seq(o_len)[1:3], o_names[1:3], elem_cl[1:3],
                              featdim[1:3], sampdim[1:3]), "\n")
                  cat(" ...")
                  cat(sprintf("\n [%i] %s: %s with %s rows and %s columns",
                              seq(o_len)[(n-2):n], o_names[(n-2):n], elem_cl[(n-2):n],
                              featdim[(n-2):n], sampdim[(n-2):n]), "\n")
              }
          })


##' @rdname QFeatures-class
##' @param x An instance of class `QFeatures`.
##' @importFrom methods callNextMethod
##' @exportMethod [
setMethod("[", c("QFeatures", "ANY", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              ## Subset the assays
              ans <- callNextMethod(x, i, j, ..., drop)
              ## Subset the AssayLinks
              ans@assayLinks <- ans@assayLinks[names(ans)]
              ## Removed lost links
              allist <- lapply(ans@assayLinks, function(al) {
                  lost <- !al@from %in% names(ans)
                  if (any(lost)) {
                      if (all(lost)) {
                          al <- AssayLink(name = al@name)
                      } else {
                          al@from <- al@from[!lost]
                          al@hits <- al@hits[!lost]
                      }
                  }
                  al
              })
              ans@assayLinks <- do.call(AssayLinks, allist)
              ## Check new object
              if (validObject(ans))
                  return(ans)
          })


##' @rdname QFeatures-class
##' @importFrom BiocGenerics dims
##' @exportMethod dims
setMethod("dims", "QFeatures",
          function(x) vapply(experiments(x), dim, integer(2)))


##' @rdname QFeatures-class
setMethod("[", c("QFeatures", "character", "ANY", "ANY"),
          function(x, i, j, k, ..., drop = TRUE) {
              if (missing(j)) j <- TRUE
              if (missing(k)) k <- TRUE
              subsetByFeature(x, i)[, j, k]
          })

##' @rdname QFeatures-class
##' 
##' @param use.names A `logical(1)` indicating whether the rownames of
##'     each assay should be propagated to the corresponding `rowData`.
##' 
setMethod("rowData", "QFeatures",
          function(x, use.names = TRUE, ...) {
              List(lapply(experiments(x), function(xx) 
                  mcols(xx, use.names = use.names, ...)))
          })

##' @rdname QFeatures-class
##'
##' @param x An instance of class `QFeatures`.
##' @param rowvars A `character()` with the names of the `rowData`
##'     variables (columns) to retain in any assay. All other
##'     variables will be dropped. In case an element in `rowvars` 
##'     isn't found in any `rowData` variable, a message is printed.
##'
##' @export
selectRowData <- function(x, rowvars) {
    stopifnot(inherits(x, "QFeatures"))
    rowvars <- as.character(rowvars)
    allvars <- unique(unlist(rowDataNames(x)))
    missingvars <- setdiff(rowvars, allvars)
    if (length(missingvars))
        message(length(missingvars), " missing/mis-typed rowvars.")
    for (i in seq_len(length(x))) {
        rd <- rowData(x[[i]])
        rowData(x[[i]]) <- rd[, colnames(rd) %in% rowvars]
    }
    x
}



##' @rdname QFeatures-class
##'     
##' @param replacement A named `list()` of same length as `i`. The 
##'     names of the list are used to select the assays to modify. Each
##'     element of the list should contain a table of class `data.frame`
##'     or `DFrame`. The names in `rowDataCols` are used to select the
##'     columns in the replacement tables and to create the new column 
##'     in the `rowData`. If instead of a table, the elements of the 
##'     list are `NULL`, then the `rowDataCols` are removed. 
##'
##' @export
setRowData <- function(object, rowDataCols, replacement) {
    stopifnot(inherits(object, "QFeatures"))
    stopifnot(is.list(replacement))
    stopifnot(all(names(replacement) %in% names(object)))
    el <- experiments(object)
    for (ii in names(replacement))
        rowData(el[[ii]])[, rowDataCols] <- replacement[[ii]][, rowDataCols]
    BiocGenerics:::replaceSlots(object,
                                ExperimentList = el,
                                check = FALSE)
}

##' @rdname QFeatures-class
##' 
##' @export
getRowData <- function(object, i, rowDataCols)  {
    ## Extract the rowData and column names from the desired assay(s)
    rdlist <- rowData(object)[i] 
    rdNames <- rowDataNames(object)[i]
    ## Check the rowDataCols argument 
    if (missing(rowDataCols)) {
        rowDataCols <- Reduce(intersect, rdNames)
        if (!length(rowDataCols)) stop("No common columns between rowData tables were found.")
    }
    ## Check that the selected columns exist in all rowData
    for (ii in seq_along(rdNames)) {
        if(!all(rowDataCols %in% rdNames[[ii]]))
            stop("Some 'rowDataCols' are not found in the rowData.")
    }
    ## Subset the rowData
    rdlist <- lapply(rdlist, function(x) x[, rowDataCols, drop = FALSE])
    rdlist <- lapply(names(rdlist), 
                     function(x) cbind(assay = x, 
                                       rowname = rownames(rdlist[[x]]),
                                       rdlist[[x]]))
    rdlist <- do.call(rbind, rdlist)
    rownames(rdlist) <- NULL
    rdlist
}


##' @rdname QFeatures-class
##'
##' @importFrom Biobase fData
##'
##' @export
rowDataNames <- function(x) {
    stopifnot(inherits(x, "MultiAssayExperiment"))
    CharacterList(lapply(experiments(x),
                         function(xx) {
                             if (inherits(xx, "SummarizedExperiment"))
                                 colnames(rowData(xx))
                             else if (inherits(xx, "eSet"))
                                 colnames(Biobase::fData(xx))
                             else NA_character_
                         }))
}


##' @rdname QFeatures-class
##'
##' @param value A character() with new name(s) for the assay(s) in `x`
##'
##' @exportMethod names<-
setReplaceMethod("names", c("QFeatures", "character"),
                 function(x, value) {
                     key_vals <- cbind(names(x), value)
                     x <-  callNextMethod(x, value)
                     names(x@assayLinks) <- value
                     for (i in seq_len(length(x))) {
                         al <- x@assayLinks[[i]]
                         al@name  <- unname(key_vals[key_vals[, 1] == al@name, 2])
                         if (!is.na(al@from))
                             al@from <- unname(key_vals[key_vals[, 1] == al@from, 2])
                         x@assayLinks[[i]] <- al
                     }
                     x
                 })


##' @rdname QFeatures-class
##'
##' @param colDataCols A `character()` that selects column(s) in the 
##'     `colData`.
##' @param rowDataCols A `character()` that selects column(s) in the
##'     `rowData`.
##' @param index The assay indicator for `SummarizedExperiment` 
##'     objects. A vector input is supported in the case that the 
##'     `SummarizedExperiment` object(s) has more than one assay 
##'     (default `1L`)
##' 
##' @importFrom MultiAssayExperiment longFormat
##' 
##' @export
longFormat <- function(object, 
                       colDataCols = NULL,
                       rowDataCols = NULL, 
                       index = 1L) {
    if (!is.null(rowDataCols)) {
        rdNames <- rowDataNames(object)
        misNames <- sapply(rdNames, 
                           function (x) any(!rowDataCols %in% x))
        ## Check that all required
        if (any(misNames))
            stop("Some 'rowDataCols' not found in assay(s): ",
                 paste0(names(misNames)[misNames], collapse = ", "))
        ## Get long format table with quantification values and colDataCols
        longDataFrame <- 
            MultiAssayExperiment::longFormat(object, colDataCols, index)
        ## Get the required rowData
        rds <- lapply(rowData(object),
                      function(rd) rd[, rowDataCols, drop = FALSE])
        rds <- do.call(rbind, rds)
        ## Merge the rowData to the long table
        cbind(longDataFrame,
              rds[longDataFrame$rowname, , drop = FALSE])
    } else {
        ## If rowDataCols is null, return the MAE longFormat output
        MultiAssayExperiment::longFormat(object, colDataCols, index)
    }
}

##' @rdname QFeatures-class
##' 
##' @param replacement A `list()` of same length as `i`. The elements
##'     of `value` must be named after `i`. Each element should 
##'     contain the replacement values to insert in the rowData. If 
##'     `length(r) > 1`, each element should be a table with the same
##'     rows as its corresponding assay and contain the `rowDataCols`
##'     variable(s). Set `replacement = NULL` to **remove** the
##'     `rowDataCols`. 
##'
##' @export
replaceRowDataCols <- function(object, 
                               replacement) {
    ## Check arguments
    stopifnot(inherits(object, "QFeatures"))
    
    el <- experiments(object)

    ## If replacement is NULL, simply remove the columns
    # Check the arguments
    stopifnot(is.list(replacement))
    stopifnot(all(names(replacement) %in% names(object)))
    ## Perform the replacement
    for (ii in names(replacement)) {
        if (nrow(rowData(el[[ii]])) != nrow(replacement[[ii]]))
            stop("'rowData' and 'replacement' don't have the same ",
                 "number of rows for assay '", ii, "'")
        rowData(el[[ii]])[, colnames(replacement[[ii]])] <- replacement[[ii]]
    }
    
    ## Since the features or columns did not change, we can safely
    ## bypass the MAE checks
    BiocGenerics:::replaceSlots(object,
                                ExperimentList = el,
                                check = FALSE)
}

##' @rdname QFeatures-class
##' 
##' @export
removeRowDataCols <- function(object, 
                              i,
                              rowDataCols) {
    ## Check arguments
    stopifnot(inherits(object, "QFeatures"))
    stopifnot(is.character(rowDataCols))
    if (missing(i)) i <- names(object)
    
    ## Remove the rowData column(s)
    el <- experiments(object)
    for (ii in i) rowData(el[[ii]])[, rowDataCols] <- NULL
    
    ## Since the features or columns did not change, we can safely
    ## bypass the MAE checks
    BiocGenerics:::replaceSlots(object,
                                ExperimentList = el,
                                check = FALSE)
}
