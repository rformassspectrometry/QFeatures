##' Quantitative MS Features
##'
##' Conceptually, a `Features` object holds a set of `matrix` (or
##' `array`) elements containing quantitative data. The number of
##' columns (samples) are always the same across the matrices, but the
##' number of rows (features) can vary. Each one of these matrices has
##' a set of feature annotations (encoded as `DataFrame` objects), that
##' have the same number of rows as the assay matrix their are
##' associated to, and an arbitrary number of columns (feature
##' variables). Such a matrix and feature annotation pair is called a
##' [FeatureSet]. In addition, a `Features` object also uses a single
##' `DataFrame` to annotate the samples (columns) represented in all
##' the matrices.
##'
##' A typical use case for such `Features` object is to represent
##' quantitative proteomics (or metabolomics) data, where different
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
##' - The features are stored as [FeatureSet] instances, and each
##'   [FeatureSet] has the same number of columns.
##'
##' - The number of rows and the names of `colData` match the number
##'   and names of the samples (columns) in the [FeatureSet]s.
##'
##' - General metadata about the object itself is provided in the
##'   `metadata` list, which can be empty.
##'
##' In the usage section, the respective arguments correspond to:
##' 
##' @param x An instance of class `Features`.
##'
##' @param object An instance of class `Features`.
##'
##' @param i A `character(1)` or `numeric(1)` for subsetting.
##'
##' @param value The replacement value.
##'
##' @section Accessors:
##'
##' The following accessors, inheritied from [SimpleList] are
##' available:
##'
##' - `length(x)`: returns the number of [FeatureSet] instances in
##'    `x`'s [FeatureSet]s. 
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
##' - `metadata(x)`, `metadata(x) <- value`: gets and sets the
##'    object's global metadata. `value` must be a `list`.
##'
##' - `x[[i]]`, `x[[i]] <- value`: gets and sets a single [FeatureSet]
##'   element, where `i` is either a `numeric(1)` or a
##'   `character(1)`. `value` must be of class `FeatureSet`.
##'
##' - `x[i]`, `x[i] <- value`: gets the list of [FeatureSet] instances
##'   defined the `numeric` or `character` `i`. Can also be used to
##'   replace the the subset of [FeatureSet] elements with a new
##'   [Features] `values` of matching samples.
##'
##' The following accessors are also available.
##'
##' - `sampleNames(object)`, `sampleNames(x) <- value`: gets and sets
##'    the samples names. `value` must be a `character` of appropriate
##'    length.
##' 
##' - `colData(x)`, `colData(x) <- value`: gets or sets the
##'    columns/samples metadata.`value` must be a `DataFrame`
##'    object. Row names of `value` match the existing column names of
##'    the assays.
##'
##' - `featureVariables(x)`: gets the list of feature variables, where
##'   each element of the list is a `character` of names for the
##'   corresponding [FeatureSet].
##'
##' - `featureVariables(object, i)`: a convenient alternative (to
##'   `featureVariables(x)[[i]]`) to get the `i`th featureData
##'   element.
##'
##' - `assays(x)`: gets a list of assays.
##'
##' - `assay(x, i)`(x, i): A convenient alternative (to
##'   `assays(x)[[i]]`) to get the `i`th (default is the largest)
##'   assay element.
##'
##' - `featureData(x)`: gets all assays' feature metadata as a `list`.
##' 
##' - `featureData(x, i)`: a convenient alternative (to
##'   `featureData(x)[[i]]`) to get the `i`th featureData element.
##' 
##' - `featureNames(x)`: gets the list of feature names, where each element of
##'    the list is a `character` of names for the corresponding assay.
##'
##' - `featureNames(object, i)`: a convenient alternative (to
##'   `featureNames(x)[[i]]`) to get the `i`th featureData element.
##' 
##'
##' @seealso The [FeatureSet] class and the [readFeatures()]
##'     constructor.
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
##' rownames(m2) <- letters[23:26]
##'
##' ## two corresponding feature metadata with appropriate row names
##' df1 <- DataFrame(Fa = 1:10, Fb = letters[1:10],
##'                  row.names = rownames(m1))
##' df2 <- DataFrame(row.names = rownames(m2))
##'
##' ## two FeatureSets, that combine the respecitve assays and feature
##' ## metadata
##' fs1 <- FeatureSet(m1, df1)
##' fs1 
##'
##' fs2 <- FeatureSet(m2, df2)
##' fs2
##'
##' ## Sample annotation (colData)
##' cd <- DataFrame(Var1 = rnorm(4),
##'                 Var2 = LETTERS[1:4],
##'                 row.names = sample_names)
##'
##' ## Putting the FeatureSets and their common sample annoation
##' ## together into a Features object
##' fts1 <- Features(fs1, fs2,
##'                  colData = cd)
##' fts1
##'
##' ## -------------------------------------------------
##' ## Creating a Features object from a data.frame (see
##' ## ?readFeatures) for details
##' ## -------------------------------------------------
##' 
##' data(hlpsms)
##' fts2 <- readFeatures(hlpsms, ecol = 1:10, name = "psms")
##' fts2
##'
##' features(fts2)
##' fts2[[1]]
##' fts2[["psms"]]
NULL

setClass("Features",
         contains = "SimpleList",
         slots = c(
             colData = "DataFrame",             
             version = "character"),
         prototype = prototype(
             elementType = "FeatureSet",
             version = "0.1"))

setMethod("show", "Features",
          function(object) {
              cat("class:", class(object), "\n")
              expt <- names(metadata(object))
              if (is.null(expt))
                  expt <- character(length(metadata(object)))
              scat("metadata(%d): %s\n", expt)
              nms <- names(object)
              if (is.null(nms))
                  nms <- character(length(object))
              scat("FeatureSets(%d): %s\n", nms)
              scat("Samples(%d): %s\n", sampleNames(object))
              scat("colData(%d): %s\n", names(colData(object)))
          })


##' @exportMethod dims
##' @rdname Features-class
setMethod("dims", "Features",
          function(x) sapply(x, dim))

##' @exportMethod dim
##' @rdname Features-class
setMethod("dim", "Features",
          function(x) dim(x[[main_assay(x)]]))

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
                     if (length(value) != ncol(object[[1]]))
                         stop("Number of sample names doesn't match")
                     rownames(object@colData) <- value
                     for (i in seq_along(object))
                         sampleNames(object[[i]]) <- value
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
          function(x, ...) lapply(x, featureVariables))

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
          function(x, i, ...) x@listData)


##' @exportMethod assays
##' @rdname Features-class
##' @param withDimnames Always `TRUE`.
setMethod("assays", "Features", function(x, ..., withDimnames = TRUE) lapply(x, assay))

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
            stop(msg, "\n'", i, "' not in names(<", class(x), ">)")
        res
    })

##' @rdname Features-class
setMethod("assay", c("Features", "missing"),
    function(x, i, ...) {
        assays <- assays(x, ...)
        if (0L == length(assays))
            stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
                 "length(assays(<", class(x), ">)) is 0'")
        assays[[main_assay(x)]]
    })



##' @exportMethod featureData
##' @rdname Features-class
setMethod("featureData", c("Features", "missing"),
          function(x, i, ...) lapply(x, featureData))

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
            stop(msg, "\n'", i, "' not in names(<", class(x), ">)")
        res
    })


##' @exportMethod featureNames
##' @rdname Features-class
setMethod("featureNames", c("Features", "missing"),
          function(x, ...) lapply(x, featureNames))

##' @rdname Features-class
setMethod("featureNames", c("Features", "numeric"),
          function(x, i, ...) rownames(featureData(x, i)))

##' @rdname Features-class
setMethod("featureNames", c("Features", "character"),
          function(x, i, ...) rownames(featureData(x, i)))

