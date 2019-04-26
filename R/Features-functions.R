## ------------------
## Private functions
## ------------------

main_assay <- function(object)
    which.max(sapply(object@featureList, nrow))

## ----------------------------
## Internal validity functions
## ----------------------------

.valid_Features_colData <- function(object) {
    if (isEmpty(object@featureList)) {
        if (nrow(object@colData) != 0)
            stop("Samples in colData but none in features")
    } else {
        n1 <- nrow(object@colData)
        n2 <- ncol(object@featureList[[1]])
        if (n1 != n2)
            stop("Number of samples in features and colData dont' match")
    }
    NULL
}

.valid_Features <- function(object) {
    .valid_Features_colData(object)    
}

setValidity("Features", .valid_Features)
