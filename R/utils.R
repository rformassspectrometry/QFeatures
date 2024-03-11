## tidy_DataFrame_columns <- function(object) {
##     for (nm in names(object)) {
##         x <- object[[nm]]
##         if (inherits(x, "List")) {
##             cat("Variable:", nm, "\n")
##             names(x) <- NULL
##             print(x)
##             cat("Do you want to:\n")
##             cat("  Leave as is      (l)\n")
##             cat("  Drop             (d)\n")
##             cat("  Summarise (mean) (s)\n")
##             k <- scan(n = 1L, what = character())
##             k <- match.arg(k, c("l", "d", "s"))
##             object[[nm]] <- switch(k,
##                                    l = object[[nm]],
##                                    d = NULL,
##                                    s = sapply(x, mean))
##         }
##     }
##     object
## }


## TODO: a programmatic version tidyDataFrame, that takes a vector of
## the form c(fvar = "d", ....) and probably one that drops all List
## columns.
##
## This begs the question whether it is necessary at all to have the
## reduced versions, as one should be able to return to the source
## FeatureSet, that contains all the feature variables.
##
## It might be useful to define some session-wide default rules?

