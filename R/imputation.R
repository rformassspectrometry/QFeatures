##' @export
imputeMethods <- function()
    c("bpca","knn", "QRILC", "MLE",
      "MinDet", "MinProb", "min", "zero",
      "mixed", "nbavg", "none")


##' @title Quantitative proteomics data imputation
##'
##' @description
##'
##' The `impute` method performs data imputation on `Features` and
##' `SummarizedExperiment` instance using a variety of methods (see
##' below).
##'
##' Users should proceed with care when imputing data and take
##' precautions to assure that the imputation produce valid results,
##' in particular with naive imputations such as replacing missing
##' values with 0.
##'
##' @details
##'
##' There are two types of mechanisms resulting in missing values in
##' LC/MSMS experiments.
##'
##' - Missing values resulting from absence of detection of a feature,
##'   despite ions being present at detectable concentrations. For
##'   example in the case of ion suppression or as a result from the
##'   stochastic, data-dependent nature of the MS acquisition
##'   method. These missing value are expected to be randomly
##'   distributed in the data and are defined as missing at random
##'   (MAR) or missing completely at random (MCAR).
##'
##' - Biologically relevant missing values resulting from the absence
##'   of the low abundance of ions (below the limit of detection of
##'   the instrument). These missing values are not expected to be
##'   randomly distributed in the data and are defined as missing not
##'   at random (MNAR).
##'
##' MNAR features should ideally be imputed with a left-censor method,
##' such as `QRILC` below. Conversely, it is recommended to use host
##' deck methods such nearest neighbours, Bayesian missing value
##' imputation or maximum likelihood methods when values are missing
##' at random.
##'
##' Currently, the following imputation methods are available. 
##'
##' - *MLE*: Maximum likelihood-based imputation method using the EM
##'   algorithm. Implemented in the `norm::imp.norm()`. function. See
##'   [norm::imp.norm()] for details and additional parameters. Note
##'   that here, `...` are passed to the [norm::em.norm()` function,
##'   rather to the actual imputation function `imp.norm`.
##'
##' - *bpca*: Bayesian missing value imputation are available, as
##'   implemented in the `pcaMethods::pca()` function. See
##'   [pcaMethods::pca()] for details and additional parameters.
##'
##' - *knn*: Nearest neighbour averaging, as implemented in the
##'   `impute::impute.knn` function. See [impute::impute.knn()]] for
##'   details and additional parameters.
##'
##' - *QRILC*: A missing data imputation method that performs the
##'   imputation of left-censored missing data using random draws from
##'   a truncated distribution with parameters estimated using
##'   quantile regression. Implemented in the
##'   `imputeLCMD::impute.QRILC`
##'   function. [imputeLCMD::impute.QRILC()] for details and
##'   additional parameters.
##'
##' - *MinDet*: Performs the imputation of left-censored missing data
##'   using a deterministic minimal value approach. Considering a
##'   expression data with *n* samples and *p* features, for each
##'   sample, the missing entries are replaced with a minimal value
##'   observed in that sample. The minimal value observed is estimated
##'   as being the q-th quantile (default `q = 0.01`) of the observed
##'   values in that sample. Implemented in the
##'   `imputeLCMD::impute.MinDet` function. See
##'   [imputeLCMD::impute.MinDet()] for details and additional
##'   parameters.
##'
##' - *MinProb*: Performs the imputation of left-censored missing data
##'   by random draws from a Gaussian distribution centred to a
##'   minimal value. Considering an expression data matrix with *n*
##'   samples and *p* features, for each sample, the mean value of the
##'   Gaussian distribution is set to a minimal observed value in that
##'   sample. The minimal value observed is estimated as being the
##'   q-th quantile (default `q = 0.01`) of the observed values in
##'   that sample. The standard deviation is estimated as the median
##'   of the feature standard deviations. Note that when estimating
##'   the standard deviation of the Gaussian distribution, only the
##'   peptides/proteins which present more than 50\% recorded values
##'   are considered. Implemented in the `imputeLCMD::impute.MinProb`
##'   function. See [imputeLCMD::impute.MinProb()] for details and
##'   additional parameters.
##'
##' - *min*: Replaces the missing values by the smallest non-missing
##'   value in the data.
##'
##' - *zero*: Replaces the missing values by 0. See also [zeroIsNA()].
##'
##' - *mixed*: A mixed imputation applying two methods (to be defined
##'   by the user as `mar` for values missing at random and `mnar` for
##'   values missing not at random, see example) on two M[C]AR/MNAR
##'   subsets of the data (as defined by the user by a `randna`
##'   logical, of length equal to nrow(object)).
##'
##' - *nbavg*: Average neighbour imputation for fractions collected
##'   along a fractionation/separation gradient, such as sub-cellular
##'   fractions. The method assumes that the fraction are ordered
##'   along the gradient and is invalid otherwise.
##'
##'   Continuous sets `NA` value at the beginning and the end of the
##'   quantitation vectors are set to the lowest observed value in the
##'   data or to a user defined value passed as argument `k`. Them,
##'   when a missing value is flanked by two non-missing neighbouring
##'   values, it is imputed by the mean of its direct neighbours. A
##'   stretch of 2 or more missing values will not be imputed. See the
##'   example below.
##'
##' - *none*: No imputation is performed and the missing values are
##'   left untouched. Implemented in case one wants to only impute
##'   value missing at random or not at random with the *mixed*
##'   method.
##'
##' The `imputeMethods()` function returns a vector with valid
##' imputation method arguments.
##'
##' @references
##'
##' Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown,
##' Trevor Hastie, Robert Tibshirani, David Botstein and Russ B.
##' Altman, Missing value estimation methods for DNA microarrays
##' Bioinformatics (2001) 17 (6): 520-525.
##'
##' Oba et al., A Bayesian missing value estimation method for gene
##' expression profile data, Bioinformatics (2003) 19 (16): 2088-2096.
##'
##' Cosmin Lazar (2015). imputeLCMD: A collection of methods for
##' left-censored missing data imputation. R package version
##' 2.0. \url{http://CRAN.R-project.org/package=imputeLCMD}.
##'
##' Lazar C, Gatto L, Ferro M, Bruley C, Burger T. Accounting for the
##' Multiple Natures of Missing Values in Label-Free Quantitative
##' Proteomics Data Sets to Compare Imputation Strategies. J Proteome
##' Res. 2016 Apr 1;15(4):1116-25. doi:
##' 10.1021/acs.jproteome.5b00981. PubMed PMID:26906401.
##'
##' @rdname impute
##'
##' @aliases impute imputeMethods impute,SummarizedExperiment-method
##'
##' @importFrom Rcpp sourceCpp
##' @useDynLib Features
##'
##' @examples
##' imputeMethods()
##'
##' data(se_na2)
##' ## table of missing values along the rows (proteins)
##' table(rowData(se_na2)$nNA)
##' ## table of missing values along the columns (samples)
##' colData(se_na2)$nNA
##'
##' ## non-random missing values
##' notna <- which(!rowData(se_na2)$randna)
##' length(notna)
##' notna
##'
##' impute(se_na2, method = "min")
##'
##' if (require("imputeLCMD")) {
##'   impute(se_na2, method = "QRILC")
##'   impute(se_na2, method = "MinDet")
##' }
##'
##' if (require("norm"))
##'   impute(se_na2, method = "MLE")
##'
##' impute(se_na2, method = "mixed",
##'        randna = rowData(se_na2)$randna,
##'        mar = "knn", mnar = "QRILC")
##' 
##' ## neighbour averaging
##' x <- se_na2[1:4, 1:6]
##' assay(x)[1, 1] <- NA ## min value
##' assay(x)[2, 3] <- NA ## average
##' assay(x)[3, 1:2] <- NA ## min value and average
##' ## 4th row: no imputation
##' assay(x)
##'
##' assay(impute(x, "nbavg"))
"impute"

##' @param object Object with missing values to be imputed.
##' @param method `character(1)` defining the imputation method. See
##'     `imputeMethods()` for available ones.
##' @param randna `logical` of length equal to `nrow(object)` defining
##'     which rows are missing at random. The other ones are
##'     considered missing not at random. Only relevant when `methods`
##'     is `mixed`.
##' @param mar Imputation method for values missing at random. See
##'     `method` above.
##' @param mnar Imputation method for values missing not at
##'     random. See `method` above.
##' @param ...
##'
##' @rdname impute
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
        if (!requireNamespace("imputeLCMD"))
            stop("Method ", method,
                 "requires the imputeLCMD package.")
    ##
    ## imputaton methods
    ##
    if (method == "knn") {
        requireNamespace("impute")        
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
        requireNamespace("impute") || stop("Package 'norm' is required.")
        s <- norm::prelim.norm(object)  ## preliminary manipulations
        th <- norm::em.norm(s, ...) ## find the MLE
        seed <- sample(.Machine$integer.max, 1)
        norm::rngseed(seed) ## set random number generator seed
        res <- norm::imp.norm(s, th, object)  ## impute missing data under the MLE
    } else if (method == "bpca"){
        requireNamespace("pcaMethods")
        nSamples <- dim(object)[2]
        .resultBPCA <- pcaMethods::pca(object,
                                       method = "bpca",
                                       nPcs = (nSamples-1),
                                       verbose = FALSE, ...)
        res <- pcaMethods::completeObs(.resultBPCA)
    } else if (method == "QRILC") {
        res <- imputeLCMD::impute.QRILC(object, ...)[[1]]
    } else if (method == "MinDet") {
        res <- imputeLCMD::impute.MinDet(object, ...)
    } else if (method == "MinProb") {
        res <- imputeLCMD::impute.MinProb(object, ...)
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

##' @export
##' @rdname impute
setMethod("impute", "SummarizedExperiment",
          function(object, method, randna, mar, mnar, ...) {
              res <- impute_matrix(assay(object), method, randna, mar, mnar, ...)
              assay(object) <- res
              object
          })
