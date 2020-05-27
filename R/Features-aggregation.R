##' @title Aggregate an assay's quantitative features
##'
##' @description
##' 
##' This function aggregates the quantitative features of an assay,
##' applying a summarisation function (`fun`) to sets of features as
##' defined by the `fcol` feature variable. The new assay's features
##' will be named based on the unique `fcol` values.
##'
##' In addition to the results of the aggregation, the newly
##' aggregated `SummarizedExperiment` assay also contains a new
##' `aggcounts` assay containing the aggregation counts matrix,
##' i.e. the number of features that were aggregated, which can be
##' accessed with the `aggcounts()` accessor.
##'
##' @param object An instance of class [Features].
##'
##' @param i The index or name of the assay which features will be
##'     aggregated the create the new assay.
##'
##' @param fcol The feature variable of assay `i` defining how to
##'     summerise the features.
##'
##' @param name A `character(1)` naming the new assay. Default is
##'     `newAssay`. Note that the function will fail if there's
##'     already an assay with `name`.
##'
##' @param fun A function used for quantitative feature
##'     aggregation. See Details for examples.
##'
##' @param ... Additional parameters passed the `fun`.
##'
##' @return A `Features` object with an additional assay.
##'
##' @details
##'
##' Aggregation is performed by a function that takes a matrix as
##' input and returns a vector of length equal to `ncol(x)`. Examples
##' thereof are
##'
##' - [MsCoreUtils::medianPolish()] to fits an additive model (two way
##'   decomposition) using Tukey's median polish_ procedure using [stats::medpolish()];
##'
##' - [MsCoreUtils::robustSummary()] to calculate a robust aggregation
##'   using [MASS::rlm()] (default);
##'
##' - [base::colMeans()] to use the mean of each column;
##'
##' - [matrixStats::colMedians()] to use the median of each column.
##'
##' - [base::colSums()] to use the sum of each column;
##'
##'
##' @section Missing quantitative values:
##' 
##' Missing quantitative values have different effect based on the
##' aggregation method employed:
##'
##' - The aggregation functions should be able to deal with missing
##'   values by either ignoring them, and propagating them. This is
##'   often done with an `na.rm` argument, that can be passed with
##'   `...`. For example, `rowSums`, `rowMeans`, `rowMedians`,
##'   ... will ignore `NA` values with `na.rm = TRUE`, as illustrated
##'   below.
##'
##' - Missing values will result in an error when using `medpolish`,
##'   unless `na.rm = TRUE` is used. Note that this option relies on
##'   implicit assumptions and/or performes an implicit imputation:
##'   when summing, the values are implicitly imputed by 0, assuming
##'   that the `NA` represent a trully absent features; when
##'   averaging, the assumption is that the `NA` represented a
##'   genuinely missing value.
##'
##' - When using robust summarisation, individual missing values are
##'   excluded prior to fitting the linear model by robust
##'   regression. To remove all values in the feature containing the
##'   missing values, use [filterNA()].
##'
##' More generally, missing values often need dedicated handling such
##' as filtering (see [filterNA()]) or imputation (see [impute()]).
##'
##' @section Missing values in the row data:
##'
##' Missing values in the row data of an assay will also impact the
##' resulting (aggregated) assay row data, as illustrated in the
##' example below. Any feature variables (a column in the row data)
##' containing `NA` values will be dropped from the aggregated row
##' data. The reasons underlying this drop are detailed in the
##' `reduceDataFrame()` manual page: only invariant aggregated rows,
##' i.e. rows resulting from the aggregation from identical variables,
##' are preserved during aggregations.
##'
##' The situation illustrated below should however only happen in rare
##' cases and should often be imputable using the value of the other
##' aggregation rows before aggregation to preserve the invariant
##' nature of that column. In cases where an `NA` is present in an
##' otherwise variant column, the column would be dropped anyway.
##' 
##' @seealso The *Features* vignette provides an extended example and
##'     the *Processing* vignette, for a complete quantitative
##'     proteomics data processing pipeline.
##' 
##' @aliases aggregateFeatures aggregateFeatures,Features-method aggcounts aggcounts,SummarizedExperiment
##'
##' @name aggregateFeatures
##'
##' @rdname Features-aggregate
##'
##' @importFrom MsCoreUtils aggregate_by_vector robustSummary countFeatures
##'
##' @examples
##'
##' ## ---------------------------------------
##' ## An example Features with PSM-level data
##' ## ---------------------------------------
##' data(feat1)
##' feat1
##'
##' ## Aggregate PSMs into peptides
##' feat1 <- aggregateFeatures(feat1, "psms", "Sequence", name = "peptides")
##' feat1
##' 
##' ## Aggregate peptides into proteins
##' feat1 <- aggregateFeatures(feat1, "peptides", "Protein", name = "proteins")
##' feat1
##'
##' assay(feat1[[1]])
##' assay(feat1[[2]])
##' aggcounts(feat1[[2]])
##' assay(feat1[[3]])
##' aggcounts(feat1[[3]])
##' 
##' ## --------------------------------------------
##' ## Aggregation with missing quantitative values
##' ## --------------------------------------------
##' data(ft_na)
##' ft_na
##'
##' assay(ft_na[[1]])
##' rowData(ft_na[[1]])
##'
##' ## By default, missing values are propagated
##' ft2 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
##' assay(ft2[[2]])
##' aggcounts(ft2[[2]])
##'
##' ## The rowData .n variable tallies number of initial rows that
##' ## were aggregated (irrespective of NAs) for all the samples. It
##' #is always >= that the sample-level aggcounts.
##' rowData(ft2[[2]])
##' 
##' ## Ignored when setting na.rm = TRUE
##' ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums, na.rm = TRUE)
##' assay(ft3[[2]])
##' aggcounts(ft3[[2]])
##' 
##' ## -----------------------------------------------
##' ## Aggregation with missing values in the row data
##' ## -----------------------------------------------
##' ## Row data results without any NAs, which includes the
##' ## Y variables
##' rowData(ft2[[2]])
##'
##' ## Missing value in the Y feature variable
##' rowData(ft_na[[1]])[1, "Y"] <- NA
##' rowData(ft_na[[1]])
##'
##' ft3 <- aggregateFeatures(ft_na, 1, fcol = "X", fun = colSums)
##' ## The Y feature variable has been dropped!
##' assay(ft3[[2]])
##' rowData(ft3[[2]])
NULL

##' @exportMethod aggregateFeatures
##' @rdname Features-aggregate
setMethod("aggregateFeatures", "Features",
          function(object, i, fcol, name = "newAssay",
                   fun = MsCoreUtils::robustSummary, ...)
              .aggregateFeatures(object, i, fcol, name, fun, ...))


.aggregateFeatures <- function(object, i, fcol, name, fun, ...) {
    if (isEmpty(object))
        return(object)
    if (name %in% names(object))
        stop("There's already an assay named '", name, "'.")
    if (missing(fcol))
        stop("'fcol' is required.")    
    if (missing(i))
        i <- main_assay(object)
    assay_i <- assay(object, i)
    rowdata_i <- rowData(object[[i]])
    if (!fcol %in% names(rowdata_i))
        stop("'fcol' not found in the assay's rowData.")
    groupBy <- rowdata_i[[fcol]]
    
    ## Message about NA values is quant/row data
    has_na <- character()
    if (anyNA(assay_i))
        has_na <- c(has_na, "quantitative")
    if (anyNA(rowdata_i, recursive = TRUE))
        has_na <- c(has_na, "row")
    if (length(has_na)) {
        msg <- paste(paste("Your", paste(has_na, collapse = " and "),
                           " data contain missing values."),
                     "Please read the relevant section(s) in the",
                     "aggregateFeatures manual page regarding the",
                     "effects of missing values on data aggregation.")
        message(paste(strwrap(msg), collapse = "\n"))
    }

    aggregated_assay <- aggregate_by_vector(assay_i, groupBy, fun, ...)
    aggcount_assay <- aggregate_by_vector(assay_i, groupBy, countFeatures)
    aggregated_rowdata <- Features::reduceDataFrame(rowdata_i, rowdata_i[[fcol]],
                                                   simplify = TRUE, drop = TRUE,
                                                   count = TRUE)

    se <- SummarizedExperiment(assays = SimpleList(assay = aggregated_assay,
                                                   aggcounts = aggcount_assay),
                               rowData = aggregated_rowdata[rownames(aggregated_assay), ])
    hits <- findMatches(rownames(aggregated_assay), groupBy)
    elementMetadata(hits)$names_to <- rowdata_i[[fcol]][hits@to]
    elementMetadata(hits)$names_from <- rownames(assay_i)[hits@to]

    assayLinks <- AssayLink(name = name,
                            from = ifelse(is.character(i), i, names(object)[i]),
                            fcol = fcol,
                            hits = hits)
    addAssay(object,
             se,
             name = name,
             assayLinks = assayLinks)

}
