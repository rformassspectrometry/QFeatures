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


.nNA <- function(x) {
    nNA <- sum(is.na(assay(x)))/prod(dim(x))
    nNA_rows <- table(rowSums(is.na(assay(x))))
    nNA_cols <- colSums(is.na(assay(x)))
    return(list(nNA = nNA, nNArows = nNA_rows, nNAcols = nNA_cols))
}

## Internal wrapper function around `.nNA` for processing multiple assays
## @param object A `QFeatures` object
## @param i One or more indices or names of the assay(s) to be processed.
.nNAi <- function(object, i) {
    if (length(object) == 1)
        return(.nNA(object[[1]]))
    res <- lapply(i,
                  function(ii) .nNA(object[[ii]]))
    .nNAasTable(object, res, i)
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

## Internal function for formating the result of nNA as a table when applied to
## multiple samples
##
## @param object A `QFeatures` object
##
## @param res A list of results obtained after applying `nNA` to multiple assays
##     of `object`
##
## @param i indices or names of the assays that were processed.
.nNAasTable <- function(object, res, i) {
    if (length(i) == 1) return(res[[1]])
    object <- object[, , i]
    names(res) <- names(object)
    ans <- vector("list", length = 3)
    names(ans) <- c("nNA", "nNArows", "nNAcols")
    ans[[1]] <- vapply(res, "[[", 1, FUN.VALUE = numeric(1))
    ans[[3]] <- t(vapply(res, "[[", 3, FUN.VALUE = numeric(3)))
    ans2 <- matrix(0,
                   ncol = 1 + nrow(colData(object)),
                   nrow = length(object))
    rownames(ans2) <- names(object)
    colnames(ans2) <- 0:nrow(colData(object))
    for (i in seq_len(length(res))) {
        x <- res[[i]]$nNArows
        ans2[i, names(x)] <- x
    }
    ans[[2]] <- ans2
    ans
}


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
##'   first element `nNA` gives the percentage of missing values; the
##'   second element `nNArows` provides a table of the number of
##'   missing values for the features (rows) of the assay(s); the
##'   third element `nNAcols` provides the number of missing values in
##'   each sample of the assay(s).
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
              for (ii in i)
                  object[[ii]] <- zeroIsNA(object[[ii]])
              object
          })

##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("QFeatures", "numeric"),
          function(object, i) zeroIsNA(object, as.integer(i)))

##' @rdname QFeatures-missing-data
setMethod("zeroIsNA", c("QFeatures", "character"),
          function(object, i) {
              for (ii in i)
                  object[[ii]] <- zeroIsNA(object[[ii]])
              object
          })


####---- infIsNA ----####


##' @exportMethod infIsNA
##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("SummarizedExperiment", "missing"),
          function(object, i) .infIsNA(object))

##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("QFeatures", "integer"),
          function(object, i) {
            for (ii in i)
              object[[ii]] <- infIsNA(object[[ii]])
            object
          })

##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("QFeatures", "numeric"),
          function(object, i) infIsNA(object, as.integer(i)))

##' @rdname QFeatures-missing-data
setMethod("infIsNA", c("QFeatures", "character"),
          function(object, i) {
            for (ii in i)
              object[[ii]] <- infIsNA(object[[ii]])
            object
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
              if (missing(i))
                  stop("'i' not provided. You must specify which assay(s) to process.")
              for (ii in i)
                  object[[ii]] <- filterNA(object[[ii]], pNA)
              object
          })
