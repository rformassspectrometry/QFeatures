##' @title Aggregate assays' quantitative features
##'
##' @description
##'
##' This function aggregates the quantitative features of one or
##' multiple assays, applying a summarisation function (`fun`) to
##' sets of features.
##' The `fcol` variable name points to a rowData column that defines
##' how to group the features during aggregate. This variable can
##' eigher be a vector (we then refer to an *aggregation by vector*)
##' or an adjacency matrix (*aggregation by matrix*).
##'
##' The rowData of the aggregated `SummarizedExperiment` assays
##' contains a `.n` variable that provides the number of parent
##' features that were aggregated.
##'
##' When aggregating with a vector, the newly aggregated
##' `SummarizedExperiment` assays also contains a new `aggcounts` assay
##' containing the aggregation counts matrix, i.e. the number of
##' features that were aggregated for each sample, which can be
##' accessed with the `aggcounts()` accessor.
##'
##' Only the rowData columns that are invariant within a group across
##' all assays will be retained in the new assays' rowData.
##'
##' @param object An instance of class [QFeatures] or [SummarizedExperiment].
##'
##' @param i A `numeric()` or `character()` indicating the index or name
##'     of one or multiple assays that will be aggregated
##'     to create one or multiple new assays.
##'
##' @param fcol A `character(1)` naming a rowdata variable (of assay
##'     `i` in case of a `QFeatures`) defining how to aggregate the
##'     features of the assays. This variable is either a `character`
##'     or a (possibly sparse) matrix. See below for details.
##'
##' @param name A `character()` naming the new assays.
##'     `name` must have the same length as i.
##'     Default is `newAssay`. Note that the function will fail if there's
##'     already an assay with `name`.
##'
##' @param fun A function used for quantitative feature
##'     aggregation. See Details for examples.
##'
##' @param ... Additional parameters passed the `fun`.
##'
##' @return A `QFeatures` object with an additional assay or a
##'  `SummarizedExperiment` object (or subclass thereof).
##'
##' @details
##'
##' Aggregation is performed by a function that takes a matrix as
##' input and returns a vector of length equal to `ncol(x)`. Examples
##' thereof are
##'
##' - [MsCoreUtils::medianPolish()] to fits an additive model (two way
##'   decomposition) using Tukey's median polish_ procedure using
##'   [stats::medpolish()];
##'
##' - [MsCoreUtils::robustSummary()] to calculate a robust aggregation
##'   using [MASS::rlm()] (default);
##'
##' - [base::colMeans()] to use the mean of each column;
##'
##' - `colMeansMat(x, MAT)` to aggregate feature by the calculating
##'    the mean of peptide intensities via an adjacency matrix. Shared
##'    peptides are re-used multiple times.
##'
##' - [matrixStats::colMedians()] to use the median of each column.
##'
##' - [base::colSums()] to use the sum of each column;
##'
##' - `colSumsMat(x, MAT)` to aggregate feature by the summing the
##'    peptide intensities for each protein via an adjacency
##'    matrix. Shared peptides are re-used multiple times.
##'
##' See [MsCoreUtils::aggregate_by_vector()] for more aggregation functions.
##'
##' @section Missing quantitative values:
##'
##' Missing quantitative values have different effects based on the
##' aggregation method employed:
##'
##' - The aggregation functions should be able to deal with missing
##'   values by either ignoring or propagating them. This is often
##'   done with an `na.rm` argument, that can be passed with
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
##' @section Using an adjacency matrix:
##'
##' When considering non-unique peptides explicitly, i.e. peptides
##' that map to multiple proteins rather than as a protein group, it
##' is convenient to encode this ambiguity explicitly using a
##' peptide-by-proteins (sparse) adjacency matrix. This matrix is
##' typically stored in the rowdata and set/retrieved with the
##' [adjacencyMatrix()] function. It can be created manually (as
##' illustrated below) or using `PSMatch::makeAdjacencyMatrix()`.
##'
##' @seealso The *QFeatures* vignette provides an extended example and
##'     the *Processing* vignette, for a complete quantitative
##'     proteomics data processing pipeline. The
##'     [MsCoreUtils::aggregate_by_vector()] manual page provides
##'     further details.
##'
##' @aliases aggregateFeatures aggregateFeatures,QFeatures-method
##'     aggcounts aggcounts,SummarizedExperiment-method
##'     adjacencyMatrix,SummarizedExperiment-method
##'     adjacencyMatrix,QFeatures-method
##'
##' @name aggregateFeatures
##'
##' @rdname QFeatures-aggregate
##'
##' @importFrom MsCoreUtils aggregate_by_vector aggregate_by_matrix robustSummary colCounts
##'
##' @importFrom methods as
##'
##' @examples
##'
##' ## ---------------------------------------
##' ## An example QFeatures with PSM-level data
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
##' ## were aggregated (irrespective of NAs) for all the samples.
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
##'
##' ## --------------------------------------------
##' ## Using a peptide-by-proteins adjacency matrix
##' ## --------------------------------------------
##'
##' ## Let's use assay peptides from object feat1 and
##' ## define that peptide SYGFNAAR maps to proteins
##' ## Prot A and B
##'
##' se <- feat1[["peptides"]]
##' rowData(se)$Protein[3] <- c("ProtA;ProtB")
##' rowData(se)
##'
##' ## This can also be defined using anadjacency matrix, manual
##' ## encoding here. See PSMatch::makeAdjacencyMatrix() for a
##' ## function that does it automatically.
##' adj <- matrix(0, nrow = 3, ncol = 2,
##'               dimnames = list(rownames(se),
##'                               c("ProtA", "ProtB")))
##' adj[1, 1] <- adj[2, 2] <- adj[3, 1:2] <- 1
##' adj
##'
##' adjacencyMatrix(se) <- adj
##' rowData(se)
##' adjacencyMatrix(se)
##'
##' ## Aggregation using the adjacency matrix
##' se2 <- aggregateFeatures(se, fcol = "adjacencyMatrix",
##'                          fun = MsCoreUtils::colMeansMat)
##'
##' ## Peptide SYGFNAAR was taken into account in both ProtA and ProtB
##' ## aggregations.
##' assay(se2)
##'
##'
##' ## Aggregation by matrix on a QFeature object works as with a
##' ## vector
##' ft <- QFeatures(list(peps = se))
##' ft <- aggregateFeatures(ft, "peps", "adjacencyMatrix", name = "protsByMat",
##'                         fun = MsCoreUtils::colMeansMat)
##' assay(ft[[2]])
##' rowData(ft[[2]])
NULL

##' @exportMethod aggregateFeatures
##' @rdname QFeatures-aggregate
setMethod("aggregateFeatures", "QFeatures",
          function(object, i, fcol, name = "newAssay",
                   fun = MsCoreUtils::robustSummary, ...) {
              if (isEmpty(object))
                  return(object)
              ## Check arguments
              if (any(present <- name %in% names(object)))
                  stop("There's already one or more assays named: '",
                       paste0(name[present], collapse = "', '"), "'.")
              i <- .normIndex(object, i)
              if (length(i) != length(name)) stop("'i' and 'name' must have same length")
              if (length(fcol) == 1) fcol <- rep(fcol, length(i))
              if (length(i) != length(fcol)) stop("'i' and 'fcol' must have same length")

              el <- experiments(object)[i]
              rowDataColsKept <- colnames(rowData(el[[i[1]]]))
              ## Aggregate each assay
              for (j in seq_along(i)) {
                  from <- i[[j]]
                  fromAssay <- el[[from]]
                  by <- fcol[[j]]
                  ## Remove already discarded columns from rowData
                  rowDataColsKept <- intersect(rowDataColsKept,
                                               colnames(rowData(fromAssay)))
                  rowData(fromAssay) <- rowData(fromAssay)[, rowDataColsKept, drop = FALSE]
                  ## Create the aggregated assay
                  el[[j]] <- aggregateFeatures(fromAssay, by, fun, ...)
                  rowDataColsKept <- colnames(rowData(el[[j]]))
                  message("\rAggregated: ", j, "/", length(i))
              }
              names(el) <- name
              for (j in name) {
                  rowDataColsKept <- intersect(rowDataColsKept,
                                               colnames(rowData(el[[j]])))
                  rowData(el[[j]]) <- rowData(el[[j]])[, rowDataColsKept, drop = FALSE]
              }
              ## Create the new QFeatures object
              for (j in seq_along(name)) {
                  object <- addAssay(object, el[[j]], name[j])
                  object <- addAssayLink(object,
                      from = i[[j]],
                      to = name[j],
                      varFrom = fcol[[j]],
                      varTo = fcol[[j]])
              }
              object
          })


##' @exportMethod aggregateFeatures
##' @rdname QFeatures-aggregate
setMethod("aggregateFeatures", "SummarizedExperiment",
          function(object, fcol, fun = MsCoreUtils::robustSummary, ...)
              .aggregateQFeatures(object, fcol, fun, ...))


.aggregateQFeatures <- function(object, fcol, fun, ...) {
    ## Copied from the PSMatch package, given that it is not available
    ## on Bioconductor yet.
    .makePeptideProteinVector <- function(m, collapse = ";") {
        stopifnot(inherits(m, "Matrix"))
        vec <- rep(NA_character_, nrow(m))
        for (i in seq_len(nrow(m)))
            vec[i] <- paste(names(which(m[i, ] != 0)), collapse = collapse)
        names(vec) <- rownames(m)
        vec
    }
    if (missing(fcol))
        stop("'fcol' is required.")
    m <- assay(object, 1)
    rd <- rowData(object)
    if (!fcol %in% names(rd))
        stop("'fcol' not found in the assay's rowData.")
    groupBy <- rd[[fcol]]

    ## Store class of assay i in case it is not a SummarizedExperiment
    ## so that the aggregated assay can be reverted to that class
    .class <- class(object)

    ## Message about NA values is quant/row data
    has_na <- character()
    if (anyNA(m))
        has_na <- c(has_na, "quantitative")
    if (anyNA(rd, recursive = TRUE))
        has_na <- c(has_na, "row")
    if (length(has_na)) {
        msg <- paste(paste("Your", paste(has_na, collapse = " and "),
                           " data contain missing values."),
                     "Please read the relevant section(s) in the",
                     "aggregateFeatures manual page regarding the",
                     "effects of missing values on data aggregation.")
        message(paste(strwrap(msg), collapse = "\n"))
    }

    if (is.vector(groupBy) & !is.list(groupBy)) { ## atomic vectors
        aggregated_assay <- aggregate_by_vector(m, groupBy, fun, ...)
        aggcount_assay <- aggregate_by_vector(m, groupBy, colCounts)
        aggregated_rowdata <- QFeatures::reduceDataFrame(rd, rd[[fcol]],
                                                         simplify = TRUE,
                                                         drop = TRUE,
                                                         count = TRUE)
        assays <- SimpleList(assay = aggregated_assay, aggcounts = aggcount_assay)
        rowdata <- aggregated_rowdata[rownames(aggregated_assay), , drop = FALSE]
    } else if (is(groupBy, "Matrix")) {
        aggregated_assay <- aggregate_by_matrix(m, groupBy, fun, ...)
        ## Remove the adjacency matrix that should be dropped anyway
        rd[[fcol]] <- NULL
        ## Temp variable for unfolding and reducing - removed later
        rd[["._vec_"]] <- .makePeptideProteinVector(groupBy)
        rd <- unfoldDataFrame(rd, "._vec_")
        aggregated_rowdata <- reduceDataFrame(rd, rd[["._vec_"]], drop = TRUE)
        aggregated_rowdata[["._vec_"]] <- NULL
        ## Count the number of peptides per protein
        .n <- apply(groupBy != 0, 2, sum)
        aggregated_rowdata[[".n"]] <- .n[rownames(aggregated_rowdata)]

        assays <- SimpleList(assay = as.matrix(aggregated_assay)) ## to discuss
        rowdata <- aggregated_rowdata[rownames(aggregated_assay), , drop = FALSE]
    } else stop("'fcol' must refer to an atomic vector or a sparse matrix.")
    se <- SummarizedExperiment(assays = assays,
                               colData = colData(object),
                               rowData = rowdata)

    ## If the input objects weren't SummarizedExperiments, then try to
    ## convert the merged assay into that class. If the conversion
    ## fails, keep the SummarizedExperiment, otherwise use the
    ## converted object (see issue #78).
    if (.class != "SummarizedExperiment")
        se <- tryCatch(as(se, .class),
                       error = function(e) se)
    return(se)
}

##' @export
##'
##' @importFrom ProtGenerics adjacencyMatrix
##'
##' @rdname QFeatures-aggregate
##'
##' @param object An instance of class `SummarizedExperiment` or
##'     `QFeatures`.
##'
##' @param adjName `character(1)` with the variable name containing
##'     the adjacency matrix. Default is `"adjacencyMatrix"`.
##'
##' @param i The index or name of the assays to extract the advaceny
##'     matrix from. All must have a rowdata variable named `adjName`.
setMethod("adjacencyMatrix", "QFeatures",
          function(object, i, adjName = "adjacencyMatrix")
              List(lapply(experiments(object)[i],
                          .adjacencyMatrix,
                          adjName = adjName)))

setMethod("adjacencyMatrix", "SummarizedExperiment",
          function(object, adjName = "adjacencyMatrix")
              .adjacencyMatrix(object, adjName))

##' @export
##'
##' @rdname QFeatures-aggregate
##'
##' @param i When adding an adjacency matrix to an assay of a
##'     `QFeatures` object, the index or name of the assay the
##'     adjacency matrix will be added to. Ignored when `x` is an
##'     `SummarizedExperiment`.
##'
##' @param value An adjacency matrix with row and column names. The
##'     matrix will be coerced to compressed, column-oriented sparse
##'     matrix (class `dgCMatrix`) as defined in the `Matrix` package,
##'     as generaled by the [sparseMatrix()] constructor.
"adjacencyMatrix<-" <- function(object, i, adjName = "adjacencyMatrix", value) {
    validAdjacencyMatrix(value)
    ## Coerse to a sparse matrix
    value <- as(value, "sparseMatrix")
    if (inherits(object, "SummarizedExperiment")) {
        if (!identical(rownames(value), rownames(object)))
            stop("Row names of the SummarizedExperiment and the adjacency matrix must match.")
        if (adjName %in% colnames(rowData(object)))
            stop("Found an existing variable ", adjName, ".")
        rowData(object)[[adjName]] <- value
        return(object)
    }
    stopifnot(inherits(object, "QFeatures"))
    if (length(i) != 1)
        stop("'i' must be of length one. Repeat the call to add a matrix to multiple assays.")
    if (is.numeric(i) && i > length(object))
        stop("Subscript is out of bounds.")
    if (is.character(i) && !(i %in% names(object)))
        stop("Assay '", i, "' not found.")
    se <- object[[i]]
    adjacencyMatrix(se, adjName = adjName) <- value
    object[[i]] <- se
    return(object)
}

.adjacencyMatrix <- function(x, adjName = "adjacencyMatrix") {
    stopifnot(adjName %in% names(rowData(x)))
    ans <- rowData(x)[[adjName]]
    if (is.null(colnames(ans)) | is.null(rownames(ans)))
        warning("The adjacency matrix should have row and column names.")
    if (!is(ans, "sparseMatrix"))
        warning("The adjacency matrix should ideally be sparse.")
    ans
}


validAdjacencyMatrix <- function(x) {
    if (is.null(colnames(x)) | is.null(rownames(x)))
        stop("The matrix must have row and column names.")
    if (any(Matrix::rowSums(x) == 0))
        stop("rowSums() == 0 detected: peptides must belong to at least one protein.")
    if (any(Matrix::colSums(x) == 0))
        stop("colSums() == 0 detected: proteins must be identified by at least one peptide.")
    TRUE
}
