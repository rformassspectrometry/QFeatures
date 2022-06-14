## ------------------
## Utility functions
## ------------------

##' @title Count Unique Features
##' 
##' @description
##' 
##' This function counts the number of unique features per sample. A 
##' grouping structure can be provided to count higher level features
##' from assays, for example counting the number of unique proteins 
##' from PSM data. 
##'
##' @param object An object of class `QFeatures`.
##'
##' @param i  A `numeric()` or `character()` vector indicating from 
##'     which assays the `rowData` should be taken.
##'
##' @param groupBy A `character(1)` indicating the variable name in 
##'     the `rowData` that contains the grouping variable, for 
##'     instance to count the unique number of peptides or proteins 
##'     expressed in each samples (column). If `groupBy` is missing, 
##'     the number of non zero elements per sample will be stored.
##'     
##' @param colDataName A `character(1)` giving the name of the new 
##'     variable in the `colData` where the number of unique features
##'     will be stored. The name cannot already exist in the 
##'     `colData`.
##'
##' @export
##'
##' @return An object of class `QFeatures`.
##' 
##' @examples 
##' data("ft_na")
##' ## Count number of (non-missing) PSMs
##' ft_na <- countUniqueFeatures(ft_na, 
##'                              i = "na", 
##'                              colDataName = "counts")
##' ft_na$counts
##' ## Count number of unique rowData feature
##' ft_na <- countUniqueFeatures(ft_na, 
##'                              i = "na", 
##'                              groupBy = "Y",
##'                              colDataName = "Y_counts")
##' ft_na$Y_counts
##' 
countUniqueFeatures <- function(object,
                                i,
                                groupBy = NULL,
                                colDataName = "count") {
    ## Check the colData does not already contain the name
    if (colDataName %in% colnames(colData(object)))
        stop("'", colDataName, "' is already present in the colData.")
    
    snames <- unlist(colnames(object)[i])
    ## Avoid that a sample is contained in 2 different assays
    if (anyDuplicated(snames))
        stop("The same sample is present in multiple assays.")
    
    ## Initialize the vector containing the feature counts
    fcounts <- vector(length = length(snames), mode = "integer")
    names(fcounts) <- snames
    
    if (is.null(groupBy)) {
        ## If no  grouping is supplied, count the non-missing
        ## features per sample
        for (ii in i) {
            fcount <- apply(assay(object[[ii]]), 2, function(x) {
                sum(!is.na(x))
            })
            fcounts[names(fcount)] <- fcount
        }
    } else {
        ## Count the number of unique entries of groupBy
        for (ii in i) {
            fcount <- apply(assay(object[[ii]]), 2, function(x) {
                length(unique(rowData(object[[ii]])[!is.na(x), groupBy]))
            })
            fcounts[names(fcount)] <- fcount
        }
    }
    
    ## Store the counts in the colData
    colData(object)[names(fcounts), colDataName] <- fcounts
    object
}



## ------------------
## Private functions
## ------------------

main_assay <- function(object)
    which.max(vapply(experiments(object),
                     nrow,
                     numeric(1)))

number_assays_in_se <- function(object)
    lengths(lapply(experiments(object), assays))

