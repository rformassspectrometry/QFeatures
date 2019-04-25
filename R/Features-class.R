##' Quantitative MS Features
##'
##' Conceptually, a `Features` object holds a set of `matrix` (or
##' `array`) elements containing quantitative data. The number of
##' columns (samples) are always the same across the matrices, but the
##' number of rows (features) can vary. Each one of these matrices has
##' a set of feature annotation (encoded as `DataFrame` objects), that
##' have the same number of rows as the assay matrix their are
##' associated to, and an arbitrary number of columns (feature
##' variables). In addition, a `Features` object also uses a single
##' `DataFrame` to annotate the samples (columns) represented in all
##' the matrices. The largest assay matrix (the one with the highest
##' number of features) is considered the main assay, from which the
##' other ones are derived by aggregating/combining several rows into
##' a single one.
##'
##' A typical use case for such `Features` object if to represent
##' quantitative proteomics or metabolomics data, where different
##' assays represent quantitation data at the PSM (the main assay),
##' peptide and protein level, and where peptide values are computed
##' from the PSM data, and the protein-level data is calculated based
##' on the peptide-level values.
##'
##' The recommended way to create `Features` objects is the use the
##' `readFeatures()` function, that creates an instance from tabular data. The
##' `Features` constructor can be used to create objects from their bare parts.
##' It iss the user's responsability to make sure that these match the class
##' validity requirements.
##'
##' @param x An instance of class `Features`.
##'
##' @param object An instance of class `Features`.
##'
##' @param i A `character(1)` or `numeric(1)` defining which assay of feature
##'     data to access.
##'
##' @param value The replacement value.
##'
##' @param ... Additional arguments.
##'
##' @param withDimnames Inherited from generic functions, but ignored here.
##'
##' @details
##'
##' A `Features` instance must comply with the following requirements:
##'
##' - Assays are matrix-like objects that all have the same number of
##'   columns. They are passed to the constructor as an `SimpleList` object
##'   names `assays`, and can be accessed with the `assays` accessor.
##'
##' - Each assay is documented by feature variable `DataFrame`. Both have the
##'   same number of rows. These feature metadata `DataFrame`s are passed to the
##'   constructor as an `SimpleList` object names `featureData` and can be
##'   accessed with the `featureData` accessor. The number of rows and the row
##'   names of each of these `DataFrame` instances must match the number of rows
##'   and rows names of the corresponding the assay. These `DataFrame` instance
##'   are mandatory but can have no columns.
##'
##' - The samples/columns of the assays are documented in a `DataFrame` named
##'   `colData` and can be accessed with the `colData` accessor. The number of
##'   rows and the row names of this `DataFrame` must match the number of
##'   columns and column names of the assays. This `DataFrame` instance is
##'   mandatory but can have no columns.
##'
##' - A `list` protiding the global metadata annotating the object as a whole,
##'   named `metadata`. It can be accessed/replaced with the `metadata`
##'   accessor/replacement method.
##'
##' @section Accessors:
##'
##' - `length(x)`: returns the number of assays in `x`.
##'
##' - `names(x)`, `names(x) <- value`: gets or sets the assays optinal
##'    names. `value` is a `character` of appropriate of length equal to
##'    `length(x)`.
##'
##' - `dims(x)`: returns the dimensions of all assays.
##'
##' - `dim(x)`: returns the dimensions of the largest assay.
##'
##' - `isEmpty(x)`: returns `TRUE` for an object without any assays, `FALSE`
##'   otherwise.
##'
##' - `metadata(x)`, `metadata(x) <- value`: gets and sets the object's global
##'    metadata. `value` must be a `list`.
##'
##' - `colData(x)`, `colData(x) <- value`: gets or sets the columns/samples
##'    metadata.`value` must be a `DataFrame` object. Row names of `value` match
##'    the existing column names of the assays.
##'
##' - `assays(x)`: gets the assays of the object. The assays are returned as a
##'   `SimpleList`.
##'
##' - `assay(x, i)`: a convenient alternative (to `assays(x)[[i]]`) to get the
##'   `i`th (default first) assay element.
##'
##' - `featureData(object)`: gets the assays' feature metadata as an instance of
##'   class `SimpleList`.
##'
##' - `featureData(object, i)`: a convenient alternative (to
##'   `featureData(x)[[i]]`) to get the `i`th featureData element.
##'
##' - `featureNames(x)`: gets the list of feature names, where each element of
##'   the list is a `character` of names for the corresponding assay.
##'
##' - `featureNames(object, i)`: a convenient alternative (to
##'   `featureNames(x)[[i]]`) to get the `i`th featureData element.
##'
##' - `featureVariables(x)`: gets the list of feature variabesl, where each
##'   element of the list is a `character` of names for the corresponding assay.
##'
##' - `featureVariables(object, i)`: a convenient alternative (to
##'   `featureVariables(x)[[i]]`) to get the `i`th featureData element.
##'
##' - `sampleNames(object)`, `sampleNames(x) <- value`: gets and sets the
##'    samples names. `value` must be a `character` of appropriate length.
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
##' ## Creating a Features object manually
##'
##' m1 <- matrix(1:40, ncol = 4)
##' m2 <- matrix(1:16, ncol = 4)
##' colnames(m1) <- colnames(m2) <- paste0("S", 1:4)
##' rownames(m1) <- letters[1:10]
##' rownames(m2) <- letters[23:26]
##'
##' df1 <- DataFrame(Fa = 1:10, Fb = letters[1:10],
##'                  row.names = letters[1:10])
##' df2 <- DataFrame(row.names = letters[23:26])
##'
##' ft1 <- Features(assays = SimpleList(assay1 = m1, assay2 = m2),
##'                 featureData = SimpleList(assay1 = df1, assay2 = df2),
##'                 colData = DataFrame(row.names = colnames(m1)),
##'                 metadata = list(paste("Generated on", Sys.Date())))
##' ft1
##'
##' ## Creating a Features object from a data.frame
##' data(hlpsms)
##' ft2 <- readFeatures(hlpsms, ecol = 1:10)
##' ft2
##'
##' ## The assay isn't named yet, so let's name it 'psms', to clarify that the
##' ## data are PSM-level quantitations.
##' names(ft2)
##' names(ft2) <- "psms"
##' ft2
NULL


setClass("Features",
    representation(
        assays = "SimpleList",
        featureData = "SimpleList",
        colData = "DataFrame",
        metadata = "list",
        version = "character"),
    prototype = prototype(version = "0.1"))

setValidity("Features", .valid_Features)

setMethod("show", "Features",
          function(object) {
              selectSome <- S4Vectors:::selectSome
              scat <- function(fmt, vals = character(), exdent = 2, ...) {
                  vals <- ifelse(nzchar(vals), vals, "''")
                  lbls <- paste(S4Vectors:::selectSome(vals), collapse = " ")
                  txt <- sprintf(fmt, length(vals), lbls)
                  cat(strwrap(txt, exdent = exdent, ...), sep="\n")
              }
              if (isEmpty(object)) .show_empty_Features(object)
              else .show_Features(object)
          })


##' @exportMethod length
##' @rdname Features-class
setMethod("length", "Features", function(x) length(x@assays))



##' @exportMethod names
##' @rdname Features-class
setMethod("names", "Features", function(x) names(x@assays))

##' @exportMethod 'names<-'
##' @rdname Features-class
setReplaceMethod("names", c("Features", "character"),
    function(x, value) {
        names(x@assays) <- value
        names(x@featureData) <- value
        x
    })

##' @exportMethod dims
##' @rdname Features-class
setMethod("dims", "Features",
          function(x) sapply(x@assays, dim))

##' @exportMethod dim
##' @rdname Features-class
setMethod("dim", "Features",
          function(x) dim(x@assays[[main_assay(x)]]))

##' @exportMethod isEmpty
##' @rdname Features-class
setMethod("isEmpty", "Features", function(x) length(x) == 0)

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

##' @exportMethod assays
##' @rdname Features-class
setMethod("assays", "Features", function(x, ...) x@assays)

##' @exportMethod assay
##' @rdname Features-class
setMethod("assay", c("Features", "numeric"),
          function(x, i, ...) {
              tryCatch({
                  assays(x, ...)[[i]]
              }, error=function(err) {
                  stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
                       "invalid subscript 'i'\n", conditionMessage(err))
              })
          })
##' @rdname Features-class
setMethod("assay", c("Features", "character"),
    function(x, i, ...) {
        msg <- paste0("'assay(<", class(x), ">, i=\"character\", ...)' ",
                      "invalid subscript 'i'")
        res <- tryCatch({
            assays(x, ...)[[i]]
        }, error=function(err) {
            stop(msg, "\n", conditionMessage(err))
        })
        if (is.null(res))
            stop(msg, "\n'", i, "' not in names(assays(<", class(x), ">))")
        res
    })

##' @rdname Features-class
setMethod("assay", c("Features", "missing"),
    function(x, i, ...) {
        assays <- assays(x, ...)
        if (0L == length(assays))
            stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
                 "length(assays(<", class(x), ">)) is 0'")
        assays[[1]]
    })


##' @exportMethod featureData
##' @rdname Features-class
setGeneric("featureData", function(x, i, ...) standardGeneric("featureData"))

##' @rdname Features-class
setMethod("featureData", "Features", function(x, ...) x@featureData)

##' @rdname Features-class
setMethod("featureData", c("Features", "numeric"),
          function(x, i, ...) {
              tryCatch({
                  featureData(x, ...)[[i]]
              }, error=function(err) {
                  stop("'featureData(<", class(x), ">, i=\"numeric\", ...)' ",
                       "invalid subscript 'i'\n", conditionMessage(err))
              })
          })

##' @rdname Features-class
setMethod("featureData", c("Features", "character"),
    function(x, i, ...) {
        msg <- paste0("'featureData(<", class(x), ">, i=\"character\", ...)' ",
                      "invalid subscript 'i'")
        res <- tryCatch({
            featureData(x, ...)[[i]]
        }, error=function(err) {
            stop(msg, "\n", conditionMessage(err))
        })
        if (is.null(res))
            stop(msg, "\n'", i, "' not in names(featureData(<", class(x), ">))")
        res
    })


##' @exportMethod featureNames
##' @rdname Features-class
setGeneric("featureNames", function(x, i, ...) standardGeneric("featureNames"))

##' @rdname Features-class
setMethod("featureNames", "Features", function(x, ...) lapply(x@featureData, rownames))

##' @rdname Features-class
setMethod("featureNames", c("Features", "numeric"),
          function(x, i, ...) rownames(featureData(x, ...)[[i]]))

##' @rdname Features-class
setMethod("featureNames", c("Features", "character"),
          function(x, i, ...) rownames(featureData(x, ...)[[i]]))



##' @exportMethod featureVariables
##' @rdname Features-class
setGeneric("featureVariables", function(x, i, ...) standardGeneric("featureVariables"))

##' @rdname Features-class
setMethod("featureVariables", "Features",
          function(x, ...) lapply(x@featureData, colnames))

##' @rdname Features-class
setMethod("featureVariables", c("Features", "numeric"),
          function(x, i, ...) colnames(featureData(x, ...)[[i]]))

##' @rdname Features-class
setMethod("featureVariables", c("Features", "character"),
          function(x, i, ...) colnames(featureData(x, ...)[[i]]))

##' @exportMethod sampleNames
##' @rdname Features-class
setMethod("sampleNames", "Features",
          function(object) rownames(object@colData))

##' @exportMethod sampleNames<-
##' @rdname Features-class
setReplaceMethod("sampleNames",c("Features", "character"),
                 function(object, value) {
                     if (isEmpty(object))
                         stop("No samples in empty object")
                     if (length(value) != ncol(assay(object)))
                         stop("Number of sample names doesn't match")
                     rownames(object@colData) <- value
                     for (i in seq_along(object@assays))
                         colnames(object@assays[[i]]) <- value
                     if (validObject(object))
                         return(object)
                 })