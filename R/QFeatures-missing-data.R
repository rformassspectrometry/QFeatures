

###---- Internal functions ----####


.zeroIsNA <- function(x) {
    sel <- assay(x) == 0
    assay(x)[sel] <- NA
    x
}

.infIsNA <- function(x) {
    sel <- is.infinite(assay(x)) 
    assay(x)[sel] <- NA
    x
}

## Internal function that compute the number and percent of missing 
## data from a SummarizedExperiment object
.nNAByAssay <- function(object) {
    x <- assay(object)
    nNA <- sum(is.na(x))
    pNA <- nNA / length(x) * 100
    DataFrame(nNA = unname(nNA), 
              pNA = unname(pNA))
}

## Internal function that compute the number and percent of missing 
## data for a given margin (features = 1 and sample = 2) from a 
## SummarizedExperiment object
.nNAByMargin <- function(object, MARGIN = 1) {
    x <- assay(object)
    nNA <- apply(is.na(x), MARGIN, sum)
    n <- ifelse(MARGIN == 1, ncol(x), nrow(x))
    pNA <- nNA / n * 100
    DataFrame(name = names(pNA), 
              nNA = unname(nNA), 
              pNA = unname(pNA))
}

## Internal function that compute the number and percent of missing 
## data for a SummarizedExperiment object
.nNA <- function(x) {
    nNA <- .nNAByAssay(x)
    nNA_rows <- .nNAByMargin(x, 1)
    nNA_cols <- .nNAByMargin(x, 2)
    list(nNA = nNA, nNArows = nNA_rows, nNAcols = nNA_cols)
}

## Internal function that compute the number and percent of missing 
## data for a QFeatures object
.nNAi <- function(object, i) {
    i <- .normIndex(object, i)
    ## Get number of missing data per assay 
    nNAassay <- do.call(rbind, lapply(i, function(ii)
        cbind(assay = ii, .nNAByAssay(object[[ii]])) ))
    ## Get number of missing data per row
    nNArow <- do.call(rbind, lapply(i, function(ii)
        cbind(assay = ii, .nNAByMargin(object[[ii]], 1)) ))
    ## Get number of missing data per column 
    nNAcol <- do.call(rbind, lapply(i, function(ii)
        cbind(assay = ii, .nNAByMargin(object[[ii]], 2)) ))
    ## Return as list
    list(nNA = nNAassay, nNArows = nNArow, nNAcols = nNAcol)
}

.row_for_filterNA <- function(x, pNA = 0L) {
    if (!is.matrix(x))
        stop(sQuote("x"), " must be a matrix.")
    if (!is.numeric(pNA))
        stop(sQuote("pNA"), " must be numeric.")
    if (length(pNA) > 1)
        stop(sQuote("pNA"), " must be of length one.")
    if (pNA > 1) pNA <- 1
    if (pNA < 0) pNA <- 0
    k <- rowSums(is.na(x)) / ncol(x)
    k <= pNA
}


###---- Documentation ----####


##' @title Managing missing data
##'
##' @description
##'
##' This manual page describes the handling of missing values in
##' [QFeatures] objects. In the following functions, if `object` is of
##' class `QFeatures`, and optional assay index or name `i` can be
##' specified to define the assay (by name of index) on which to
##' operate.
##'
##' The following functions are currently available:
##'
##' - `zeroIsNA(object, i)` replaces all 0 in `object` by `NA`. This
##'    is often necessary when third-party software assume that
##'    features that weren't quantified should be assigned an
##'    intensity of 0.
##'
##' - `infIsNA(object, i)` replaces all infinite values in `object` by 
##'    `NA`. This is necessary when third-party software divide 
##'    expression data by zero values, for instance during custom 
##'    normalization.
##'
##' - `nNA(object, i)` return a list of missing value summaries. The
##'   first element `nNA` gives a `DataFrame` with the number and the
##'   percentage of missing values for the whole assay; the second 
##'   element `nNArows` provides a `DataFrame` of the number and the 
##'   percentage of missing values for the features (rows) of the 
##'   assay(s); the third element `nNAcols` provides the number and 
##'   the percentage of missing values in each sample of the assay(s).
##'   When `object` has class `QFeatures` and additional column with 
##'   the assays is provided in each element's `DataFrame`.
##'
##' - `filterNA(object, pNA, i)` removes features (rows) that contain
##'   `pNA` percentage or more missing values.
##'
##' See the *Processing* vignette for examples.
##'
##' @param object An object of class `QFeatures` or `SummarizedExperiment`.
##'
##' @param pNA `numeric(1)` providing the maximim percentage of
##'     missing values per feature (row) that is acceptable. Feature
##'     with higher percentages are removed. If 0 (default), features
##'     that contain any number of `NA` values are dropped.
##'
##' @param i One or more indices or names of the assay(s) to be processed.
##'
##' @return An instance of the same class as `object`.
##'
##' @aliases zeroIsNA zeroIsNA,SummarizedExperiment,missing-method zeroIsNA,QFeatures,missing-method zeroIsNA,QFeatures,numeric-method zeroIsNA,QFeatures,integer-method zeroIsNA,QFeatures,character-method
##'
##' @aliases infIsNA infIsNA,SummarizedExperiment,missing-method infIsNA,QFeatures,missing-method infIsNA,QFeatures,numeric-method infIsNA,QFeatures,integer-method infIsNA,QFeatures,character-method
##'
##' @aliases nNA nNA,SummarizedExperiment,missing-method nNA,QFeatures,missing-method nNA,QFeatures,numeric-method nNA,QFeatures,integer-method nNA,QFeatures,character-method
##'
##' @aliases filterNA filterNA,SummarizedExperiment-method filterNA,QFeatures-method
##'
##' @name missing-data
##'
##' @rdname QFeatures-missing-data
##'
##' @seealso The `impute()` for `QFeautres` instances.
##'
##' @examples
##' se_na2
##'
##' ## Summary if missing values
##' nNA(ft_na, 1)
##'
##' ## Remove rows with missing values
##' assay(filterNA(ft_na, i = 1))
##'
##' ## Replace NAs by zero and back
##' ft_na <- impute(ft_na, i = 1, method = "zero")
##' assay(ft_na)
##' ft_na <- zeroIsNA(ft_na, 1)
##' assay(ft_na)
NULL


####---- zeroIsNA ----####


##' @exportMethod zeroIsNA
##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("SummarizedExperiment", "missing"),
          function(object, i) .zeroIsNA(object))

##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("QFeatures", "integer"),
          function(object, i) {
              el <- experiments(object)
              for (ii in i)
                  el[[ii]] <- zeroIsNA(el[[ii]])
              replaceAssay(object, el)
          })

##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("QFeatures", "numeric"),
          function(object, i) zeroIsNA(object, as.integer(i)))

##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("QFeatures", "character"),
          function(object, i) {
              if (any(! i %in% names(object)))
                  stop("subscript contains invalid names")
              zeroIsNA(object, which(names(object) %in% i))
          })


####---- infIsNA ----####


##' @exportMethod infIsNA
##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("SummarizedExperiment", "missing"),
          function(object, i) .infIsNA(object))

##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("QFeatures", "integer"),
          function(object, i) {
              el <- experiments(object)
              for (ii in i)
                  el[[ii]] <- infIsNA(el[[ii]])
              replaceAssay(object, el)
          })

##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("QFeatures", "numeric"),
          function(object, i) infIsNA(object, as.integer(i)))

##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("QFeatures", "character"),
          function(object, i) {
              if (any(! i %in% names(object)))
                  stop("subscript contains invalid names")
              infIsNA(object, which(names(object) %in% i))
          })


####---- nNA ----####


##' @exportMethod nNA
##' @rdname QFeatures-missing-data
setMethod("nNA", c("SummarizedExperiment", "missing"),
          function(object, i) .nNA(object))

##' @rdname QFeatures-missing-data
setMethod("nNA", c("QFeatures", "integer"),
          function(object, i) .nNAi(object, i))

##' @rdname QFeatures-missing-data
setMethod("nNA", c("QFeatures", "numeric"),
          function(object, i) .nNAi(object, as.integer(i)))

##' @rdname QFeatures-missing-data
setMethod("nNA", c("QFeatures", "character"),
          function(object, i) .nNAi(object, i) )


####---- filterNA ----####


##' @exportMethod filterNA
##' @rdname QFeatures-missing-data
setMethod("filterNA", "SummarizedExperiment",
          function(object, pNA = 0) {
              k <- .row_for_filterNA(assay(object), pNA)
              object[k, ]
          })

##' @rdname QFeatures-missing-data
setMethod("filterNA", "QFeatures",
          function(object, pNA = 0, i) {
              if (!is.character(i)) i <- names(object)[i]
              sel <- lapply(i, function(ii) {
                  .row_for_filterNA(assay(object[[ii]]), pNA)
              })
              names(sel) <- i 
              object[sel, ]
          })
