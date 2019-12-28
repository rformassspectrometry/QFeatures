imputeMethods <- function()
    c("bpca","knn", "QRILC", "MLE",
      "MinDet", "MinProb", "min", "zero",
      "mixed", "nbavg", "none")


impute_matrix <- function(object,
                          method,
                          randna,
                          mar,
                          mnar,
                          ...) {
    if (!anyNA(object)) return(object)
    if (missing(method))
        stop("Please specify an imputation method. ",
             "See '?impute' for details.")
    method <- match.arg(method,
                        choices = imputeMethods(),
                        several.ok = FALSE)
    res <- object
    if (method %in% c("CRILC", "MinDet", "MinProb"))
        if (!require("imputeLCMD"))
            stop("Method ", method,
                 "requires the imputeLCMD package.")
    ##
    ## imputaton methods
    ##
    if (method == "knn") {
        imp_res <- impute::impute.knn(object, ...)
        res <- imp_res$data
        if (!is.null(imp_res$rng.state)) {
            assign(".Random.seed", imp_res$rng.state, envir = .GlobalEnv)
        } else {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    } else if (method == "nbavg") {
        message("Assuming values are ordered.")
        impargs <- pairlist(...)
        if (is.null(impargs$k)) k <- min(object, na.rm = TRUE)
        else k <- impargs$k
        res <- imp_neighbour_avg(object, k = k)
    } else if (method == "MLE") {
        require("norm") || stop("Package 'norm' is required.")
        s <- norm::prelim.norm(object)  ## preliminary manipulations
        th <- norm::em.norm(s, ...) ## find the MLE
        seed <- sample(.Machine$integer.max, 1)
        norm::rngseed(seed) ## set random number generator seed
        res <- norm::imp.norm(s, th, x)  ## impute missing data under the MLE
    } else if (method == "bpca"){
        nSamples <- dim(object)[2]
        .resultBPCA <- pca(object, method = "bpca",
                           nPcs = (nSamples-1), verbose = FALSE, ...)
        res <- completeObs(.resultBPCA)
    } else if (method == "QRILC") {
        res <- imputeLCMD::impute.QRILC(exprs(object), ...)[[1]]
    } else if (method == "MinDet") {
        res <- imputeLCMD::impute.MinDet(exprs(object), ...)
    } else if (method == "MinProb") {
        res <- imputeLCMD::impute.MinProb(exprs(object), ...)
    } else if (method == "min") {
        val <- min(object, na.rm = TRUE)
        res[is.na(res)] <- val
    } else if (method == "mixed") {
        if (missing(randna))
            stop("Mixed imputation requires 'randna' argument. See ?impute.",
                 call. = FALSE)
        stopifnot(is.logical(randna))
        if (missing(mar))
            stop("Mixed imputation requires 'mar' argument. See ?impute.",
                 call. = FALSE)
        if (missing(mnar))
            stop("Mixed imputation requires 'mnar' argument. See ?impute.",
                 call. = FALSE)
        if (length(randna) != nrow(res))
            stop("Number of proteins and length of randna must be equal.",
                 call. = FALSE)
        res[randna, ] <- impute_matrix(object[randna, ], mar, ...)
        res[!randna, ] <- impute_matrix(object[!randna, ], mnar, ...)
    } else if (method == "zero") {
        res[is.na(res)] <- 0
    }
    ## else method == "none" -- do nothing
    res
}
