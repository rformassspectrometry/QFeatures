##' @title hyperLOPIT PSM-level expression data
##'
##' @description
##'
##' A `data.frame` with PSM-level quantitation data by Christoforou *et al.*
##' (2016). This is the first replicate of a spatial proteomics dataset from a
##' hyperLOPIT experimental design on Mouse E14TG2a embryonic stem
##' cells. Normalised intensities for proteins for TMT 10-plex labelled
##' fractions are available for 3 replicates acquired in MS3 mode using an
##' Orbitrap Fusion mass-spectrometer.
##'
##' The variable names are
##'
##' - X126, X127C, X127N, X128C, X128N, X129C, X129N, X130C, X130N and
##'   X131: the 10 TMT tags used to quantify the peptides
##'   along the density gradient.
##'
##' - Sequence: the peptide sequence.
##'
##' - ProteinDescriptions: the description of the protein this peptide was
##'   associated to.
##'
##' - NbProteins: the number of proteins in the protein group.
##'
##' - ProteinGroupAccessions: the main protein accession number in the protein
##'   group.
##'
##' - Modifications: post-translational modifications identified in the peptide.
##'
##' - qValue: the PSM identification q-value.
##'
##' - PEP: the PSM posterior error probability.
##'
##' - IonScore: the Mascot ion identification score.
##'
##' - NbMissedCleavages: the number of missed cleavages in the peptide.
##'
##' - IsolationInterference: the calculated precursor ion isolation interference.
##'
##' - IonInjectTimems: the ions injection time in milli-seconds.
##'
##' - Intensity: the precursor ion intensity.
##'
##' - Charge: the peptide charge.
##'
##' - mzDa: the peptide mass to charge ratio, in Daltons.
##'
##' - MHDa: the peptide mass, in Daltons.
##'
##' - DeltaMassPPM: the difference in measure and calculated mass, in parts per
##'   millions.
##'
##' - RTmin: the peptide retention time, in minutes.
##'
##' - markers: localisation for well known sub-cellular markers. QFeatures of
##'   unknown location are encode as `"unknown"`.
##'
##' For further details, install the `pRolocdata` package and see
##' `?hyperLOPIT2015`.
##'
##' @source
##'
##' The `pRolocdata` package: \url{http://bioconductor.org/packages/pRolocdata/}
##'
##' @references
##'
##' *A draft map of the mouse pluripotent stem cell spatial proteome*
##' Christoforou A, Mulvey CM, Breckels LM, Geladaki A, Hurrell T, Hayward PC,
##' Naake T, Gatto L, Viner R, Martinez Arias A, Lilley KS. Nat Commun. 2016 Jan
##' 12;7:8992. doi: 10.1038/ncomms9992. PubMed PMID: 26754106; PubMed Central
##' PMCID: PMC4729960.
##'
##' @seealso
##'
##' See [QFeatures] to import this data using the [readQFeatures()] function.
##'
##' @md
"hlpsms"


##' Feature example data
##'
##' `feat1` is a small test `QFeatures` object for testing and
##' demonstration. `feat2` is used to demonstrate assay joins. `ft_na`
##' is a tiny test set that contains missing values used to
##' demonstrate and test the impact of missing values on data
##' processing. `se_na2` is an `SummarizedExperiment` with missing
##' values of mixed origin.
##'
##' @aliases ft_na se_na2 feat2
##'
"feat1"

##' Example `QFeatures` object after processing
##' 
##' `feat3` is a small `QFeatures` object that contains 7 assays: 
##' `psms1`, `psms2`, `psmsall`, `peptides`, `proteins`, 
##' `normpeptides`, `normproteins`. The dataset contains example data
##' that could be obtained after running a simple processing pipeline.
##' You can use it to get your hands on manipulating `AssayLinks` 
##' since all 3 general cases are present:
##' * One parent to one child `AssayLink`: the relationship can either be 
##'   one row to one row (e.g. "peptides" to "normpeptides") or 
##'   multiple rows to one row (e.g. "peptides" to "proteins").
##' * One parent to multiple children `AssayLink`: for instance "peptides"
##'   to "normpeptides" and "proteins".
##' * Multiple parents to one child `AssayLink`: links the rows between 
##'   multiple assays to a single assays where some rows in different 
##'   parent assays may point to the same row in the child assay. E.g.
##'   "psms1" and "psms2" to "psmsall"
##'
##' @source
##'
##' `feat3` was built from `feat1`. The source code is available in
##' [`inst/scripts/test_data.R`](https://github.com/rformassspectrometry/QFeatures/blob/master/inst/scripts/test_data.R)
##' 
##' @seealso
##'
##' See `?feat1` for other example/test data sets. 
##' 
##' @examples 
##' 
##' data("feat3")
##' plot(feat3)
##' 
##' @md
"feat3"