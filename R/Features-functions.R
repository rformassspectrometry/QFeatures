.main_assay <- function(object)
    which.max(sapply(object@assays, nrow))

.show_empty_Features <- function(object)
    cat("Empty", class(object), "object\n")

.show_Features <- function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")

    ## metadata()
    expt <- names(metadata(object))
    if (is.null(expt))
        expt <- character(length(metadata(object)))
    scat("metadata(%d): %s\n", expt)

    ## assays()
    nms <- assayNames(object)
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
    scat("assays(%d): %s\n", nms)
    scat("Features(%d): %s\n", featureNames(object))
    scat("Samples(%d): %s\n", sampleNames(object))

    ## colData()
    scat("Feature variables(%d): %s\n", names(colData(object)))
}


.valid_Features_assay_dims <- function(object) {
    ## checking assays dimensions
    all_dims <- sapply(object@assays, function(m) dim(m)[1:2])
    if (any(is.na(all_dims)))
        return(wmsg("all assays must be matrix-like objects ",
                    "with 2 (or more?) dimensions"))
    if (!all(all_dims[2L, ] == all_dims[2L, 1L]))
        stop("all assays must have the same number of cols")

    NULL
}

.valid_Features_colData <- function(object) {
    ## checking colData dimensions
    p1 <- ncol(object@assays[[1]])
    p2 <- nrow(object@colData)
    if (p1 != p2)
        stop("Number of samples and columns metadata don't match")

    NULL
}

.valid_Features_fData <- function(object) {
    ## checking feature data dimension
    n1 <- max(sapply(object@assays, nrow))
    n2 <- nrow(object@fData)
    if (n1 != n2)
        stop("Number of features and feature metadata don't match")

    NULL
}

.valid_Features_rownames <- function(object) {
    ## checking row names
    rn1 <- rownames(object@assays[[.main_assay(object)]])
    rn2 <- rownames(object@fData)
    if (!identical(rn1, rn2))
        stop("Feature names don't match in assay and feature metadata")

    NULL
}


.valid_Features <- function(object) {
    if (isEmpty(object))
        return(NULL)
    .valid_Features_assay_dims(object)
    .valid_Features_colData(object)
    .valid_Features_fData(object)
    .valid_Features_rownames(object)
}
