
## ----------------------------
## Internal validity functions
## ----------------------------

##' @importFrom methods slot
.valid_QFeatures_indices <- function(object) {
    if (!isEmpty(object) && !identical(names(object), names(object@assayLinks)))
        stop("Assay links names are wrong.")
    NULL
}

.unique_row_names <- function(object) {
    dup_row_names <- vapply(experiments(object),
                            function(x) anyDuplicated(rownames(x)),
                            numeric(1))
    if (any(dup_row_names != 0))
        stop("Assay(s) ", paste(which(dup_row_names != 0), collapse = ", "),
             " has/have duplicated row names.")
    NULL
}

## The function checks that features that are linked from or to in a
## Hits object are present in the corresponding linked assay. 
## @param hits A `Hits` object to check the links
## @param linkedAssay The assay linked from or to (an object that 
##     inherits from `SummarizedExperiment`)
## @param direction Either "from" or "to"
.checkLinksInHits <- function(hits, linkedAssay, direction = "from") {
    ## Get the feature names in the `Hits` object
    linkedFNames <- mcols(hits)[, paste0("names_", direction)]
    if (is.null(linkedAssay) & length(linkedFNames) > 0) {
        ## If there is no linkedAssay, there cannot be links in Hits
        stop("@hits contains links that point ", direction, 
             " a missing assay")
    } else if (!all(linkedFNames %in% rownames(linkedAssay))) {
        ## All links in Hits must point to or from existing features
        ## in the linkedAssay
        stop("@hits contains links that point ", direction, 
             " missing features")
    }
    NULL
}

## This function checks whether a given AssayLink in the QFeatures
## object is valid, if not it return an informative error. 
## 
## @param object A QFeatures object 
## @param i A `numeric(1)`, `character(1)` or `logical(1)` providing 
##     the index of the assayLink to check
.validAssayLink <- function(object, i) {
    al <- assayLink(object, i)
    
    ## An AssayLink object is valid if the name(s) of the parent 
    ## assay(s) for all assays are either NA (= root node) or 
    ## contained in the assay names of the QFeatures object
    alFrom <- al@from
    if (!all(is.na(alFrom) | alFrom %in% names(object)))
        stop("@from not valid")
    
    ## An AssayLink object is valid if the links encoded in `@hits` do
    ## link from and to existing features in the `QFeatures` object.
    if (inherits(al@hits, "List")) { ## If the AssayLink contains a HitsList object 
        ## Get the parent assays. Note that this sophisticated line of
        ## code allows that if any of the elements in `al@from` is NA,
        ## it returns a NULL element required for later checks
        parents <- lapply(al@from, function(ii) object[[ii]])
        ## Check links from parents
        mapply(function(hits, linkedAssay) {
            .checkLinksInHits(hits, linkedAssay, direction = "from")
        }, hits = al@hits, linkedAssay = parents)
        ## Check links to self
        lapply(al@hits, .checkLinksInHits, 
               linkedAssay = object[[al@name]], direction = "to")
    } else { ## If the AssayLink contains a single Hits object
        ## Check links from parent
        .checkLinksInHits(al@hits, object[[al@from]], direction = "from")
        ## Check links to self
        .checkLinksInHits(al@hits, object[[al@name]], direction = "to")
    }
    NULL
}

## This function checks whether the AssayLinks in the QFeatures object
## are valid, if not it return an informative error. 
## 
## @param object a QFeatures object 
## 
.validAssayLinks <- function(object) {
    ## An AssayLinks object is valid if the names of the node for all 
    ## assays are contained in the assay names of the QFeatures object.
    ## The order matters!
    assayn <- names(object)
    aln <- vapply(object@assayLinks, "slot", "name", 
                  FUN.VALUE = character(1))
    if (!all(aln == assayn))
        stop("@names not valid")
    
    ## An AssayLinks object is valid if all its AssayLink objects are
    ## valid
    for (i in names(object)) {
        .validAssayLink(object, i)
    }
    NULL
}


.valid_QFeatures <- function(object) {
    .valid_QFeatures_indices(object)
    .unique_row_names(object)
    .validAssayLinks(object)
}

## ----------------------------
## Export the validity function
## ----------------------------

setValidity("QFeatures", .valid_QFeatures)
