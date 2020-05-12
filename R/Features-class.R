##' @title Quantitative MS Features
##'
##' @description
##'
##' Conceptually, a `Features` object holds a set of *assays*, each
##' composed of a `matrix` (or `array`) containing quantitative data
##' and row annotations (meta-data).  The number and the names of the
##' columns (samples) must always be the same across the assays, but
##' the number and the names of the rows (features) can vary. The
##' assays are typically defined as `SummarizedExperiment` objects. In
##' addition, a `Features` object also uses a single `DataFrame` to
##' annotate the samples (columns) represented in all the matrices.
##'
##' The `Features` class extends the
##' [MultiAssayExperiment::MultiAssayExperiment] and inherits all
##' the functionality of the
##' [MultiAssayExperiment::MultiAssayExperiment] class.
##'
##' A typical use case for such `Features` object is to represent
##' quantitative proteomics (or metabolomics) data, where different
##' assays represent quantitation data at the PSM (the main assay),
##' peptide and protein level, and where peptide values are computed
##' from the PSM data, and the protein-level data is calculated based
##' on the peptide-level values. The largest assay (the one with the
##' highest number of features, PSMs in the example above) is
##' considered the main assay.
##'
##' The recommended way to create `Features` objects is the use the
##' `readFeatures()` function, that creates an instance from tabular
##' data. The `Features` constructor can be used to create objects
##' from their bare parts.  It is the user's responsability to make
##' sure that these match the class validity requirements.
##'
##' @section Constructors:
##'
##' - `Features(..., assayLinks)` allows the manual construction of
##'   objects. It is the user's responsability to make sure these
##'   comply. The arguments in `...` are those documented in
##'   [MultiAssayExperiment::MultiAssayExperiment()]. For details
##'   about `assayLinks`, see [AssayLinks]. An example is shown below.
##'
##' - The [readFeatures()] function constructs a `Features` object
##'   from text-based spreadsheet or a `data.frame` used to generate
##'   an assay. See the function manual page for details and an
##'   example.
##'
##' @section Accessors:
##'
##' - The `Features` class extends the
##'   [MultiAssayExperiment::MultiAssayExperiment] class and inherits
##'   all its accessors and replacement methods.
##'
##' - The `rowDataNames` accessor returns a list with the `rowData`
##'   variable names.
##'
##' @section Adding assays:
##'
##' - The [aggregateFeatures()] function creates a new assay by
##'   aggregating features of an existing assay.
##'
##' - `addAssay(x, y, name, assayLinks)`: Adds a new assay (or
##'   list of assays) `y` to the `Features` instance `x`. `name`
##'   is a `character(1)` naming the single assay (default is
##'   `"newAssay"`), and is ignored if `y` is a list of
##'   assays. `assayLinks` is an optional [AssayLinks].
##'
##' @section Subsetting:
##'
##' - Features object can be subset using the `x[i, j, k, drop =
##'   TRUE]` paradigm. See the argument descriptions for details.
##'
##' - The [subsetByFeature()] function can be used to subset a
##'   `Features` object using one or multiple feature names that will
##'   be matched across different assays, taking the aggregation
##'   relation between assays.
##'
##' - The `selectRowData(x, rowvars)` function can be used to
##'   select a limited number of `rowData` columns of interest named
##'   in `rowvars` in the `x` instance of class `Features`.
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
##' @seealso
##'
##' - The [readFeatures()] constructor and the [aggregateFeatures()]
##'   function. The *Features* vignette provides an extended example.
##'
##' - The [Features-filtering] manual page demonstrates how to filter
##'   features based on their rowData.
##'
##' - The [missing-data] manual page to manage missing values in
##'   `Features` objects.
##'
##' - The [Features-processing] and [aggregateFeatures()] manual pages
##'   and *Processing* vignette describe common quantitative data
##'   processing methods using in quantitative proteomics.
##'
##' @import MultiAssayExperiment ProtGenerics
##'
##' @name Features
##'
##' @aliases Features Features-class class:Features addAssay dims,Features-method show,Features-method [,Features,ANY,ANY,ANY-method [,Features,character,ANY,ANY-method
##'
##' @aliases rowDataNames selectRowData
##'
##' @rdname Features-class
##'
##' @exportClass Features
##'
##' @author Laurent Gatto
##'
##' @examples
##' ## ------------------------
##' ## An empty Features object
##' ## ------------------------
##'
##' Features()
##'
##' ## -----------------------------------
##' ## Creating a Features object manually
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
##' fts1 <- Features(el, colData = cd)
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
##' ## Keep only the Fa variable
##' selectRowData(fts1, rowvars = "Fa")
##'
##' ## -----------------------------------
##' ## See ?readFeatures to create a
##' ## Features object from a data.frame
##' ## or spreadsheet.
##' ## -----------------------------------
NULL


## ----------------------------------
## Features Class ChangeLog
##
## Version 0.1:
##  - Contains MatchedAssayExperiment
## Version 0.2:
##  - Contains MultiAssayExperiment (see issue 46)

setClass("Features",
         contains = "MultiAssayExperiment",
         slots = c(version = "character",
                   assayLinks = "AssayLinks"),
         prototype = prototype(
             version = "0.2"))


##' @rdname Features-class
##' @param object An instance of class `Features`.
##' @exportMethod show
setMethod("show", "Features",
          function(object) {
              if (isEmpty(object)) {
                  cat(sprintf("A empty instance of class %s", class(object)), "\n")
                  return(NULL)
              }
              cat(sprintf("A instance of class %s", class(object)), "containing")
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
              cat(sprintf("\n [%i] %s: %s with %s rows and %s columns",
                          seq(o_len), o_names, elem_cl, featdim, sampdim), "\n")
          })



##' @rdname Features-class
##' @param x An instance of class `Features`.
##' @importFrom methods callNextMethod
##' @exportMethod [
setMethod("[", c("Features", "ANY", "ANY", "ANY"),
          function(x, i, j, ..., drop = TRUE) {
              ## Subset the assays
              ans <- callNextMethod(x, i, j, ..., drop)
              ## Subset the AssayLinks
              ans@assayLinks <- ans@assayLinks[names(ans)]
              ## Removed lost links
              allist <- lapply(ans@assayLinks, function(al) {
                  if (!al@from %in% names(ans)) al <- AssayLink(name = al@name)
                  al
              })
              ans@assayLinks <- do.call(AssayLinks, allist)
              ## Check new object
              if (validObject(ans))
                  return(ans)
          })


##' @rdname Features-class
##' @importFrom BiocGenerics dims
##' @exportMethod dims
setMethod("dims", "Features",
          function(x) sapply(experiments(x), dim))


##' @rdname Features-class
setMethod("[", c("Features", "character", "ANY", "ANY"),
          function(x, i, j, k, ..., drop = TRUE) {
              if (missing(j)) j <- TRUE
              if (missing(k)) k <- TRUE
              subsetByFeature(x, i)[, j, k]
          })

##' @rdname Features-class
##'
##' @param x An instance of class `Features`.
##' @param rowvars A `character()` with the names of the `rowData`
##'     variables (columns) to retain in any assay. All other
##'     variables will be dropped. In case an element in `rowvars`
##'     isn't found in any `rowData` variable, a message is printed.
##'
##' @export
selectRowData <- function(x, rowvars) {
    stopifnot(inherits(x, "Features"))
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


##' @rdname Features-class
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


##' @rdname Features-class
##' 
##' @param value A character() with new name(s) for the assay(s) in `x`  
##' 
##' @exportMethod names<-
setReplaceMethod("names", c("Features", "character"),
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
