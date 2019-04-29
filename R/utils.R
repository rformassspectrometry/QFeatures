tidy_DataFrame_columns <- function(object) {
    for (nm in names(object)) {
        x <- object[[nm]]
        if (inherits(x, "List")) {
            cat("Variable:", nm, "\n")
            names(x) <- NULL
            print(x)
            cat("Do you want to:\n")
            cat("  Leave as is      (l)\n")
            cat("  Drop             (d)\n")
            cat("  Summarise (mean) (s)\n")
            k <- scan(n = 1L, what = character())
            k <- match.arg(k, c("l", "d", "s"))
            object[[nm]] <- switch(k,
                                   l = object[[nm]],
                                   d = NULL,
                                   s = sapply(x, mean))
        }
    }
    invisible(object)
}


##' @rdname tidyFeatureData
##' @exportMethod tidyFeatureData
setMethod("tidyFeatureData", c("FeatureSet", "missing"),
          function(object, i, ...) {
              df <- featureData(object)
              df <- tidy_DataFrame_columns(df)
              object@featureData <- df
              if (validObject(object))
                  object
          })



##' Manually review feature variables
##'
##' This function offers a user to either keep, drop or summarise a
##' feature variable to amend a [FeatureSet]'s or a [Features]'
##' `featureData` slot. This is necessary when proceeding with
##' multiple calls to `combineFeatures` to avoid overly intricated
##' 'lists of lists' variables.
##'
##' @param object An instance of class [FeatureSet] or [Features].
##' 
##' @param i For [Features] objects only, defines which `featureData`
##'     to review.
##'
##' @param ... Ignored.
##'
##' @return An updated instance of the same class as `object.
##'
##' @md
##' @exportMethod tidyFeatureData
##' @rdname tidyFeatureData
##' 
##' @aliases tidyFeatureData 
##' 
setMethod("tidyFeatureData", c("Features", "character"),
          function(object, i, ...) {
              x <- object[[i]]
              x <- tidyFeatureData(x)
              object@listData[[i]] <- x
              if (validObject(object))
                  object
          })

##' @exportMethod tidyFeatureData
##' @rdname tidyFeatureData
setMethod("tidyFeatureData", c("Features", "numeric"),
          function(object, i, ...) {
              x <- object[[i]]
              x <- tidyFeatureData(x)
              object@listData[[i]] <- x
              if (validObject(object))
                  object
          })


