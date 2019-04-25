##' @title Quantitative MS Features
##'
##' @aliases Features Features-class class:Features
##'
##' @name Features
##'
##' @description
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
##' `[readFeatures()]` function, that creates an instance from tabular data. The
##' `Features` constructor can be used to create objects from their bare parts.
##' It iss the user's responsability to make sure that these match the class
##' validity requirements.
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
##' @section Accessors and
##'
##' In the following section `x` and `object` are `Features` objects.
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
##' - `sampleNames(object)`: gets the samples names.
##'
##' @rdname Features-class
##'
##' @import S4Vectors
##' @importFrom SummarizedExperiment assays assay colData colData<-
##' @importFrom Biobase sampleNames
##'
##' @md
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
setMethod("length", "Features", function(x) length(x@assays))

##' @exportMethod names
setMethod("names", "Features", function(x) names(x@assays))
##' @exportMethod 'names<-'
setReplaceMethod("names", c("Features", "character"),
    function(x, value) {
        names(x@assays) <- value
        names(x@featureData) <- value
        x
    })

##' @exportMethod dims
setMethod("dims", "Features",
          function(x) sapply(x@assays, dim))

##' @exportMethod dim
setMethod("dim", "Features",
          function(x) dim(x@assays[[main_assay(x)]]))

##' @exportMethod isEmpty
setMethod("isEmpty", "Features", function(x) length(x) == 0)

##' @exportMethod metadata
setMethod("metadata", "Features", function(x, ...) x@metadata)

##' @exportMethod metadata<-
setReplaceMethod("metadata", c("Features", "list"),
    function(x, value) {
        if (!length(value))
            names(value) <- NULL
        x@metadata <- value
        x
    })

##' @exportMethod colData
setMethod("colData", "Features", function(x) x@colData)

##' @exportMethod colData<-
setReplaceMethod("colData", c("Features", "DataFrame"),
    function(x, value) {
        x@colData <- value
        if (validObject(x)) x
    })

##' @exportMethod assays
setMethod("assays", "Features", function(x, ...) x@assays)

##' @exportMethod assay
setMethod("assay", c("Features", "numeric"),
          function(x, i, ...) {
              tryCatch({
                  assays(x, ...)[[i]]
              }, error=function(err) {
                  stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
                       "invalid subscript 'i'\n", conditionMessage(err))
              })
          })

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

setMethod("assay", c("Features", "missing"),
    function(x, i, ...) {
        assays <- assays(x, ...)
        if (0L == length(assays))
            stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
                 "length(assays(<", class(x), ">)) is 0'")
        assays[[1]]
    })


##' @exportMethod featureData
setGeneric("featureData", function(x, i, ...) standardGeneric("featureData"))

setMethod("featureData", "Features", function(x, ...) x@featureData)

setMethod("featureData", c("Features", "numeric"),
          function(x, i, ...) {
              tryCatch({
                  featureData(x, ...)[[i]]
              }, error=function(err) {
                  stop("'featureData(<", class(x), ">, i=\"numeric\", ...)' ",
                       "invalid subscript 'i'\n", conditionMessage(err))
              })
          })

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
setGeneric("featureNames", function(x, i, ...) standardGeneric("featureNames"))

setMethod("featureNames", "Features", function(x, ...) lapply(x@featureData, rownames))

setMethod("featureNames", c("Features", "numeric"),
          function(x, i, ...) rownames(featureData(x, ...)[[i]]))

setMethod("featureNames", c("Features", "character"),
          function(x, i, ...) rownames(featureData(x, ...)[[i]]))



##' @exportMethod featureVariables
setGeneric("featureVariables", function(x, i, ...) standardGeneric("featureVariables"))

setMethod("featureVariables", "Features",
          function(x, ...) lapply(x@featureData, colnames))

setMethod("featureVariables", c("Features", "numeric"),
          function(x, i, ...) colnames(featureData(x, ...)[[i]]))

setMethod("featureVariables", c("Features", "character"),
          function(x, i, ...) colnames(featureData(x, ...)[[i]]))

##' @exportMethod sampleNames
setMethod("sampleNames", "Features",
          function(object) rownames(object@colData))
