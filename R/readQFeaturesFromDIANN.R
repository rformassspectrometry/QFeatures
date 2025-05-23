##' @title  Read DIA-NN output as a QFeatures objects
##'
##' @description
##'
##' This function takes the `Report.tsv` output files from DIA-NN and
##' converts them into a multi-set `QFeatures` object. It is a wrapper
##' around [readQFeatures()] with default parameters set to match
##' DIA-NN label-free and plexDIA report files: default `runCol` is
##' `"File.Name"` and default quantCols` is `"Ms1.Area"`.
##'
##' @inheritParams readQFeatures
##'
##' @param extractedData A `data.frame` or any object that can be
##'     coerced to a `data.frame` that contains the data from the
##'     `*_ms1_extracted.tsv` file generated by DIA-NN. This argument
##'     is optional and is currently only applicable for mTRAQ
##'     multiplexed experiments where DIA-NN was run using the
##'     `plexdia` module (see references).
##'
##' @param multiplexing A `character(1)` indicating the type of
##'     multiplexing used in the experiment. One of `"none"` (default,
##'     for label-free experiments) or `"mTRAQ"` (for plexDIA
##'     experiments).
##'
##' @param ... Further arguments passed to [readQFeatures()].
##'
##' @return An instance of class `QFeatures`. The quantiative data of
##'     each acquisition run is stored in a separate set as a
##'     `SummarizedExperiment` object.
##'
##' @seealso
##'
##' - The `QFeatures` (see [QFeatures()]) class to read about how to
##'   manipulate the resulting `QFeatures` object.
##'
##' - The [readQFeatures()] function which this one depends on.
##'
##' @references
##' Derks, Jason, Andrew Leduc, Georg Wallmann, R. Gray Huffman,
##' Matthew Willetts, Saad Khan, Harrison Specht, Markus Ralser,
##' Vadim Demichev, and Nikolai Slavov. 2022. "Increasing the
##' Throughput of Sensitive Proteomics by plexDIA." Nature
##' Biotechnology, July.
##' [Link to article](http://dx.doi.org/10.1038/s41587-022-01389-w)
##'
##' @author Laurent Gatto, Christophe Vanderaa
##'
##' @importFrom tidyr pivot_wider
##' @importFrom tidyselect all_of
##'
##' @export
##'
##' @examples
##'
##' x <- read.delim(MsDataHub::benchmarkingDIA.tsv())
##' x[["File.Name"]] <- x[["Run"]]
##'
##' #################################
##' ## Label-free multi-set case
##'
##' ## using default arguments
##' readQFeaturesFromDIANN(x)
##'
##' ## use the precursor identifier as assay rownames
##' readQFeaturesFromDIANN(x, fnames = "Precursor.Id") |>
##'     rownames()
##'
##' ## with a colData (and default arguments)
##' cd <- data.frame(sampleInfo = LETTERS[1:24],
##'                  quantCols = "Ms1.Area",
##'                  runCol = unique(x[["File.Name"]]))
##' readQFeaturesFromDIANN(x, colData = cd)
##'
##' #################################
##' ## mTRAQ multi-set case
##'
##' x2 <- read.delim(MsDataHub::Report.Derks2022.plexDIA.tsv())
##' x2[["File.Name"]] <- x2[["Run"]]
##' readQFeaturesFromDIANN(x2, multiplexing = "mTRAQ")
readQFeaturesFromDIANN <- function(assayData,
                                   colData = NULL,
                                   quantCols = "Ms1.Area",
                                   runCol = "File.Name",
                                   multiplexing = c("none", "mTRAQ"),
                                   extractedData = NULL,
                                   ecol = NULL,
                                   verbose = TRUE,
                                   ...) {
    multiplexing <- match.arg(multiplexing, several.ok = FALSE)
    if (multiplexing == "mTRAQ") {
        if (verbose) message("Pivoting quantiative data.")
        assayData <- .formatMtraqReportData(assayData, colData,
                                      quantCols, runCol)
        quantCols <- assayData[[2]]
        assayData <- assayData[[1]]
    } ## else is none
    ans <- readQFeatures(assayData, colData = colData,
                         quantCols = quantCols,
                         runCol = runCol,
                         ecol = ecol,
                         ...)
    if (!is.null(extractedData)) {
        ans <- .addDiannExtractedData(ans, extractedData)
    }
    setQFeaturesType(ans, "bulk")
}

## (Only for mTRAQ multiplexing!) Internal function that extracts the
## mTRAQlabels from the peptide sequence, removes the mTRAQ annotation
## from the precursor ID, identifies constant columns within precursor
## and puts the quantification data for different mTRAQ labels in
## separate columns (wide format).
.formatMtraqReportData <- function(assayData, colData, quantCols, runCol) {
    assayData$Label <-
        sub("^.*[Q-](\\d).*$", "\\1", assayData$Modified.Sequence)
    assayData$Precursor.Id <-
        gsub("\\(mTRAQ.*?\\)", "(mTRAQ)", assayData$Precursor.Id)
    ## .checkLabelsInColData(colData, assayData)
    idCols <- .findPrecursorVariables(assayData,
                                      precursorId = "Precursor.Id",
                                      runCol = runCol)
    ans <- pivot_wider(
        assayData, id_cols = all_of(idCols),
        names_from = "Label",
        values_from = all_of(quantCols)
    )
    list(assayData = ans,
         quantCols = setdiff(colnames(ans), colnames(assayData)))
}

## Internal function that identifies which variables in the report
## data are constant within each precursor within each run.
.findPrecursorVariables <- function(assayData, precursorId, runCol) {
    precIds <- paste0(assayData[[precursorId]], assayData[[runCol]])
    nUniqueIds <- length(unique(precIds))
    nLevels <- sapply(colnames(assayData), function(x) {
        nrow(unique(assayData[, c(precursorId, runCol, x)]))
    })
    names(nLevels)[nLevels == nUniqueIds]
}

## Internal function that adds the extractedData to a QFeatures
## object. The functions first converts the extractedData to a
## SummarizedExperiment objects and subsets the SCE for the set of
## shared samples. The added assay is automatically linked (using
## AssayLinks) to the assayData.
## Developer's note: the function assumes that DIA-NN creates sample
## names in the extracted data by appending the labels to the run
## names
## @param object A QFeatures object containing DIA-NN report data, as
##     generated by readQFeatures.
## @param extractedData A data.frame or any object that can be coerced
##     to a data.frame that contains the data from the
##     `*_ms1_extracted.tsv` file generated by DIA-NN.
.addDiannExtractedData <- function(object, extractedData) {
    quantColPattern <- paste0(unique(object$Label), "$", collapse = "|")
    quantCols <- grep(quantColPattern, colnames(extractedData))
    extractedData <- readSummarizedExperiment(
        extractedData, quantCols = quantCols, fnames = "Precursor.Id"
    )
    extractedData <- .keepSharedSamples(extractedData, object)
    object <- addAssay(object, extractedData, name = "Ms1Extracted")
    addAssayLink(
        object,
        from = grep("Ms1Extracted", names(object), invert = TRUE),
        to = "Ms1Extracted",
        varFrom = rep("Precursor.Id", length(names(object)) - 1),
        varTo = "Precursor.Id"
    )
}

## Internal functions that subsets the extractedData to keep only
.keepSharedSamples <- function(extractedData, object) {
    cnames <- unique(unlist(colnames(object)))
    if (any(mis <- !cnames %in% colnames(extractedData)))
        stop("Some columns present in assayData are not found in ",
             "extracted data", paste0(cnames[mis], collapse = ", "),
             "\nAre you sure the two tables were generated from ",
             "the same experiment?")
    extractedData[, cnames]
}
