## Desired interface
##
## Question: Subset my data `x` by row (protein) `P12345`
##
## Desired output: I want the assays that contain features for that
## protein. These include:
##
## - the `proteins` assay, and specifically the feature names `P12345`.
## - the `peptides` assay, and specifically the peptides associated to
##   protein `P12345`.
## - the `psms` assay, and specifically the PSMs associated to the
##   peptides associated to protein `P12345`.
##
## Desired interface: x["P12345", ]
##
## This requires to
##
## 1. Find the assay that contains feature P12345 (this is the
##    `proteins` assay).
##
## 2. Recursively find the parent assays (these are `peptides` and
##    `psms`).
##
## 3. Find the features in these parent assay that are
##    linked/associated to the protein of interest
##


find_assay_with_feature_name <- function(x, i) {
    rnms <- rownames(x)
    which(sapply(rnms, function(x) i %in% x))
}
