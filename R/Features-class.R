##' Quantitative MS Features
##'
##' Conceptually, a `Features` object holds a set of `matrix` (or
##' `array`) elements containing quantitative data. The number of
##' columns (samples) are always the same across the matrices, but the
##' number of rows (features) can vary. Each one of these matrices has
##' a set of feature annotation (encoded as `DataFrame` objects), that
##' have the same number of rows as the assay matrix their are
##' associated to, and an arbitrary number of columns (feature
##' variables). Such a matrix and feature annotation pair is called a
##' [FeatureSet] and all these [FeatureSet] instances are available as
##' a single [FeatureList] object. In addition, a `Features` object
##' also uses a single `DataFrame` to annotate the samples (columns)
##' represented in all the matrices.
##'
##' A typical use case for such `Features` object is to represent
##' quantitative proteomics or metabolomics data, where different
##' assays represent quantitation data at the PSM (the main assay),
##' peptide and protein level, and where peptide values are computed
##' from the PSM data, and the protein-level data is calculated based
##' on the peptide-level values.
##'
##' The largest [FeatureSet] (the one with the highest number of
##' features) is considered the main assay, from which the other ones
##' are derived by aggregating/combining several rows into a single
##' one.
##'
##' The recommended way to create `Features` objects is the use the
##' `readFeatures()` function, that creates an instance from tabular
##' data. The `Features` constructor can be used to create objects
##' from their bare parts.  It is the user's responsability to make
##' sure that these match the class validity requirements.
##'
##' A `Features` instance must comply with the following requirements:
##'
##' - The features are stored as a [FeatureList], where each
##'   [FeatureSet] has the same number of columns.
##'
##' - The number of rows and the names of `colData` match the number
##'   and names of the samples (columns) in the [FeatureList]
##'   elements.
##'
##' - General metadata about the object itself is provided in the
##'   `metadata` list.
##'
##' @param x An instance of class `Features`.
##'
##' @param object An instance of class `Features`.
##'
##' @param i A `character(1)` or `numeric(1)` for subsetting.
##'
##' @param j A `character(1)` or `numeric(1)` for subsetting
##'     (currently ignored in `Features`).
##'
##' @param value The replacement value.
##'
##' @param ... Additional arguments.
##'
##' @param withDimnames Inherited from generic functions, but ignored here.
##'
##' @section Accessors:
##'
##' The following accessors are available.
##'
##' - `length(x)`: returns the number of [FeatureSet] instances in
##'    `x`'s [FeatureList].
##'
##' - `names(x)`, `names(x) <- value`: gets or sets the assays
##'    optional names of the [FeatureSet] instances. `value` is a
##'    `character` of length equal to `length(x)`.
##'
##' - `dims(x)`: returns the dimensions of all [FeatureSet]'s assays.
##'
##' - `dim(x)`: returns the dimensions of the largest assay.
##'
##' - `isEmpty(x)`: returns `TRUE` for an object without any
##'   [FeatureSet]s, `FALSE` otherwise.
##'
##' - `metadata(x)`, `metadata(x) <- value`: gets and sets the object's global
##'    metadata. `value` must be a `list`.
##'
##' - `sampleNames(object)`, `sampleNames(x) <- value`: gets and sets the
##'    samples names. `value` must be a `character` of appropriate length.
##' 
##' - `colData(x)`, `colData(x) <- value`: gets or sets the columns/samples
##'    metadata.`value` must be a `DataFrame` object. Row names of `value` match
##'    the existing column names of the assays.
##'
##' - `featureVariables(x)`: gets the list of feature variabels, where
##'   each element of the list is a `character` of names for the
##'   corresponding [FeatureSet].
##'
##' - `featureVariables(object, i)`: a convenient alternative (to
##'   `featureVariables(x)[[i]]`) to get the `i`th featureData
##'   element.
##'
##' @name Features
##'
##' @rdname Features-class
##'
##' @aliases Features Features-class class:Features
##'
##' @md
##'
##' @import S4Vectors
##' @importFrom SummarizedExperiment assays assay colData colData<-
##' @importFrom Biobase sampleNames sampleNames<-
##' @importFrom BiocGenerics dims
##'
##' @exportClass Features
##'
##' @author Laurent Gatto
##'
##' @examples
##'
##' ## An empty Features object
##' 
##' Features()
##'
##' ## Creating a Features object manually
##'
##' m1 <- matrix(1:40, ncol = 4)
##' m2 <- matrix(1:16, ncol = 4)
##' sample_names <- paste0("S", 1:4)
##' colnames(m1) <- colnames(m2) <- sample_names
##' rownames(m1) <- letters[1:10]
##' rownames(m2) <- letters[23:26]
##'
##' df1 <- DataFrame(Fa = 1:10, Fb = letters[1:10],
##'                  row.names = letters[1:10])
##' df2 <- DataFrame(row.names = letters[23:26])
##'
##' fs1 <- FeatureSet(m1, df1)
##' fs1 
##'
##' fs2 <- FeatureSet(m2, df2)
##' fs2
##'
##' fl <- FeatureList(fs1, fs2)
##' names(fl) <- paste0("Fs", 1:2)
##' fl
##'
##' fts1 <- Features(featureList = fl,
##'                  colData = DataFrame(Var = rnorm(4),
##'                                      row.names = sample_names))
##' fts1
##'
##' ## Creating a Features object from a data.frame
##' data(hlpsms)
##' fts2 <- readFeatures(hlpsms, ecol = 1:10, name = "psms")
##' fts2
##'
##' features(fts2)
##' features(fts2)[[1]]
NULL

setClass("Features",
         slots = c(featureList = "FeatureList",
                   colData = "DataFrame",
                   metadata = "list",
                   version = "character"),
         prototype = prototype(version = "0.1"))


setMethod("show", "Features",
          function(object) {
              cat("class:", class(object), "\n")
              expt <- names(metadata(object))
              if (is.null(expt))
                  expt <- character(length(metadata(object)))
              scat("metadata(%d): %s\n", expt)
              nms <- names(object)
              if (is.null(nms))
                  nms <- character(length(object@featureList))
              scat("FeatureSets(%d): %s\n", nms)
              scat("Samples(%d): %s\n", sampleNames(object))
              scat("colData(%d): %s\n", names(colData(object)))
          })

##' @exportMethod length
##' @rdname Features-class
setMethod("length", "Features", function(x) length(x@featureList))

##' @exportMethod names
##' @rdname Features-class
setMethod("names", "Features", function(x) names(x@featureList))

##' @exportMethod 'names<-'
##' @rdname Features-class
setReplaceMethod("names", c("Features", "character"),
    function(x, value) {
        names(x@featureList) <- value
        x
    })

##' @exportMethod dims
##' @rdname Features-class
setMethod("dims", "Features",
          function(x) sapply(x@featureList, dim))

##' @exportMethod dim
##' @rdname Features-class
setMethod("dim", "Features",
          function(x) dim(x@featureList[[main_assay(x)]]))

##' @exportMethod isEmpty
##' @rdname Features-class
setMethod("isEmpty", "Features", function(x) isEmpty(x@featureList))

##' @exportMethod metadata
##' @rdname Features-class
setMethod("metadata", "Features", function(x, ...) x@metadata)

##' @exportMethod metadata<-
##' @rdname Features-class
setReplaceMethod("metadata", c("Features", "list"),
    function(x, value) {
        if (!length(value))
            names(value) <- NULL
        x@metadata <- value
        x
    })

##' @exportMethod sampleNames
##' @rdname Features-class
setMethod("sampleNames", "Features",
          function(object) rownames(object@colData))

##' @exportMethod sampleNames<-
##' @rdname Features-class
setReplaceMethod("sampleNames", c("Features", "character"),
                 function(object, value) {
                     if (isEmpty(object))
                         stop("No samples in empty object")
                     if (length(value) != ncol(object@featureList[[1]]))
                         stop("Number of sample names doesn't match")
                     rownames(object@colData) <- value
                     for (i in seq_along(object@featureList))
                         sampleNames(object@featureList[[i]]) <- value
                     if (validObject(object))
                         return(object)
                 })

##' @exportMethod colData
##' @rdname Features-class
setMethod("colData", "Features", function(x) x@colData)

##' @exportMethod colData<-
##' @rdname Features-class
setReplaceMethod("colData", c("Features", "DataFrame"),
    function(x, value) {
        x@colData <- value
        if (validObject(x)) x
    })

##' @exportMethod featureVariables
##' @rdname Features-class
setMethod("featureVariables", c("Features", "missing"),
          function(x, ...) lapply(x@featureList, featureVariables))

##' @rdname Features-class
setMethod("featureVariables", c("Features", "numeric"),
          function(x, i, ...) featureVariables(x)[[i]])

##' @rdname Features-class
setMethod("featureVariables", c("Features", "character"),
          function(x, i, ...) {
              if (!i %in% names(x))
                  stop(i, "' not in names(<", class(x), ">)")
              featureVariables(x)[[i]]
          })

##' @exportMethod features
##' @rdname Features-class
setMethod("features", c("Features", "missing"),
          function(x, i, ...) x@featureList)

##' @exportMethod features
##' @rdname Features-class
setMethod("features", c("Features", "numeric"),
          function(x, i, ...) {
              tryCatch({
                  features(x, ...)[[i]]
              }, error = function(err) {
                  stop("'features(<", class(x), ">, i=\"numeric\", ...)' ",
                       "invalid subscript 'i'\n", conditionMessage(err))
              })
          })

##' @exportMethod features
##' @rdname Features-class
setMethod("features", c("Features", "character"),
    function(x, i, ...) {
        msg <- paste0("'features(<", class(x), ">, i=\"character\", ...)' ",
                      "invalid subscript 'i'")
        res <- tryCatch({
            features(x, ...)[[i]]
        }, error = function(err) {
            stop(msg, "\n", conditionMessage(err))
        })
        if (is.null(res))
            stop(msg, "\n'", i, "' not in names(features(<", class(x), ">))")
        res
    })

##' @exportMethod [
##' @rdname Features-class
setMethod("[", c("Features", "ANY", "missing"),
          function(x, i, j, ..., drop = FALSE) {
              if (missing(i))
                  return(x)
              if (!length(i))
                  return(Features())
              if (is.character(i) && (!i %in% names(x)))
                  stop(paste0("<", class(x), ">[i,] index out of bounds: %s"))
              if (is.numeric(i) && (i > nrow(colData(x))))
                  stop(paste0("<", class(x), ">[i,] index out of bounds: %s"))
              x@featureList <- x@featureList[i]
              if (validObject(x))
                  x
          })

## ---------------------------------------------------------------
## - check from here ----------------------------------------------
## ---------------------------------------------------------------

## ##' @exportMethod assays
## ##' @rdname Features-class
## setMethod("assays", "Features", function(x, ...) x@assays)

## ##' @exportMethod assay
## ##' @rdname Features-class
## setMethod("assay", c("Features", "numeric"),
##           function(x, i, ...) {
##               tryCatch({
##                   assays(x, ...)[[i]]
##               }, error=function(err) {
##                   stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
##                        "invalid subscript 'i'\n", conditionMessage(err))
##               })
##           })
## ##' @rdname Features-class
## setMethod("assay", c("Features", "character"),
##     function(x, i, ...) {
##         msg <- paste0("'assay(<", class(x), ">, i=\"character\", ...)' ",
##                       "invalid subscript 'i'")
##         res <- tryCatch({
##             assays(x, ...)[[i]]
##         }, error=function(err) {
##             stop(msg, "\n", conditionMessage(err))
##         })
##         if (is.null(res))
##             stop(msg, "\n'", i, "' not in names(assays(<", class(x), ">))")
##         res
##     })

## ##' @rdname Features-class
## setMethod("assay", c("Features", "missing"),
##     function(x, i, ...) {
##         assays <- assays(x, ...)
##         if (0L == length(assays))
##             stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
##                  "length(assays(<", class(x), ">)) is 0'")
##         assays[[1]]
##     })


## ##' @exportMethod featureData
## ##' @rdname Features-class
## setMethod("featureData", "Features", function(x, ...) x@featureData)

## ##' @rdname Features-class
## setMethod("featureData", c("Features", "numeric"),
##           function(x, i, ...) {
##               tryCatch({
##                   featureData(x, ...)[[i]]
##               }, error=function(err) {
##                   stop("'featureData(<", class(x), ">, i=\"numeric\", ...)' ",
##                        "invalid subscript 'i'\n", conditionMessage(err))
##               })
##           })

## ##' @rdname Features-class
## setMethod("featureData", c("Features", "character"),
##     function(x, i, ...) {
##         msg <- paste0("'featureData(<", class(x), ">, i=\"character\", ...)' ",
##                       "invalid subscript 'i'")
##         res <- tryCatch({
##             featureData(x, ...)[[i]]
##         }, error=function(err) {
##             stop(msg, "\n", conditionMessage(err))
##         })
##         if (is.null(res))
##             stop(msg, "\n'", i, "' not in names(featureData(<", class(x), ">))")
##         res
##     })


## ##' @exportMethod featureNames
## ##' @rdname Features-class
## setMethod("featureNames", "Features", function(x, ...) lapply(x@featureData, rownames))

## ##' @rdname Features-class
## setMethod("featureNames", c("Features", "numeric"),
##           function(x, i, ...) rownames(featureData(x, ...)[[i]]))

## ##' @rdname Features-class
## setMethod("featureNames", c("Features", "character"),
##           function(x, i, ...) rownames(featureData(x, ...)[[i]]))


## - `assays(x)`: gets the assays of the object. The assays are returned as a
##   `SimpleList`.
##
##- `assay(x, i)`: a convenient alternative (to `assays(x)[[i]]`) to get the
##   `i`th (default first) assay element.
##
##- `featureData(object)`: gets the assays' feature metadata as an instance of
##   class `SimpleList`.
##
##- `featureData(object, i)`: a convenient alternative (to
##   `featureData(x)[[i]]`) to get the `i`th featureData element.
##
##- `featureNames(x)`: gets the list of feature names, where each element of
##   the list is a `character` of names for the corresponding assay.
##
##- `featureNames(object, i)`: a convenient alternative (to
##   `featureNames(x)[[i]]`) to get the `i`th featureData element.
