.main_assay <- function(object)
    which.max(sapply(object@assays, nrow))

.valid_Features_assay_dims <- function(object) {
    l1 <- length(object@assays)
    l2 <- length(object@featureData)
    if (!identical(l1, l2))
        stop("Different number of assays and feature data")        
    all_dims <- sapply(object@assays, function(m) dim(m)[1:2])
    if (any(is.na(all_dims)))
        return(wmsg("all assays must be matrix-like objects ",
                    "with at least 2 dimensions"))
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

.valid_Features_featureData <- function(object) {
    ## checking feature data dimension
    n1 <- sapply(object@assays, nrow)
    n2 <- sapply(object@featureData, nrow)
    if (!all(n1 == n2))
        stop("Number of features in assays and feature data don't match")

    NULL
}

.valid_Features_rownames <- function(object) {
    ## checking row names
    rn1 <- lapply(object@assays, rownames)
    rn2 <- lapply(object@featureData, rownames)        
    if (!identical(rn1, rn2))
        stop("Feature names in assay and feature data don't match")
    NULL
}

.valid_Features_names <- function(object) {
    nms1 <- names(object@assays)
    nms2 <- names(object@featureData)
    if (!identical(nms1, nms2))
        stop("Names in assay and feature data don't match")
}

.valid_Features <- function(object) {
    if (isEmpty(object))
        return(NULL)
    .valid_Features_assay_dims(object)
    .valid_Features_colData(object)
    .valid_Features_featureData(object)
    .valid_Features_rownames(object)
    .valid_Features_names(object)
}


.show_empty_Features <- function(object)
    cat("Empty", class(object), "object\n")

.show_Features <- function(object) {
    scat <- function(fmt, vals=character(), exdent=2, ...) {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse = " ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep = "\n")
    }
    
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")

    ## metadata()
    expt <- names(metadata(object))
    if (is.null(expt))
        expt <- character(length(metadata(object)))
    scat("metadata(%d): %s\n", expt)

    ## assays()
    nms <- names(object)
    if (is.null(nms))
        nms <- character(length(object@assays))
    scat("assays(%d): %s\n", nms)
    scat("Features(%d): %s\n", featureNames(object))
    scat("Samples(%d): %s\n", sampleNames(object))

    ## colData()
    scat("Feature variables(%d): %s\n", featureVariables(object))
}


.featureVariables <- function(object, assay = NULL) {
    stopifnot(inherits(object, "Features"))
    if (isEmpty(object))
        return(NA_character_)
    if (is.null(assay))
        return(unlist(unique(sapply(object@featureData, names))))
    if (is.character(assay)) {
        assay <- assay[1]
        stopifnot(assay %in% names(object))
    }
    names(object@featureData[[assay]])
}
