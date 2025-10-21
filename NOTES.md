# One or two levels in QFeatures (closed)

One or two assay levels could be considered in QFeatures:

- one level: each SE contains only a single assay, and when an SE is
  processed (log-transformed, normalised, …) in a way that its
  dimensions stay the same, a new SE is created and added to the
  QFeatures object.

- two level: SEs can contain multiple assays, and when an SE is
  processed (log-transformed, normalised, …) in a way that its
  dimensions stay the same, a new assay is added to that SE.

This
[question](https://stat.ethz.ch/pipermail/bioc-devel/2020-January/016096.html)
on the bioc-devel list ask for advice on SE processing, and whether a
new SE or new assay in the original SE should be preferred. While the
letter is arguably more elegant, and is also used in
SingleAssayExperiment pipelines, it doesn’t seem to be the case when
using SummarizedExperiments.

As for features (or [MultiAssayExperiments in
general](https://github.com/waldronlab/MultiAssayExperiment/issues/266)),
the two-level approach isn’t readily available out-of-the-box, and would
require additional developments:

- Every function that operates on an SE of a QFeatures object would need
  to allow the user to specify which assay to use (and/or by default use
  the latest one).

- The `show,QFeatures` method would need to display the number/names of
  the assays in each SE to make these two levels explicit.

Despite the elegant of the two-level option, it seems that the
additional development isn’t warranted at this time.

The [`updateAssay`
function](https://github.com/rformassspectrometry/QFeatures/issues/37)
was originally intended for the two-level approach, i.e. to add an assay
to an SE. This is not considered anymore (for now, at least).

There is one exception though. When aggregating features with
[`aggregateFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md),
a second assay is added, named `aggcounts` that counts the number of
features that were aggregate for each sample and each low-level
features.

# How to add new assays (closed)

1.  Through aggregation with `aggregateFeatures`.

2.  Processing an SE.

This can/could be done explicitly with `addAssay`

    addAssay(cptac, logTransform(cptac[["peptides"]]), name = "peptides_log")
    addAssay(cptac, logTransform(cptac[[1]]), name = "peptides_log")

or implicitly

    logTransform(cptac, "peptides", name = "peptides_log")
    logTransform(cptac, 1, name = "peptides_log")

3.  Joining SEs (for example multiple TMT batches) (TODO)

&nbsp;

    joinAssays(QFeatures, c("pep_batch1", "pep_batch2", "pep_batch3"), name = "peptides")
    joinAssays(QFeatures, c(1, 2, 3), name = "peptides")

See below.

# QFeatures API

### Processing functions

- A processing function that acts on a Feature’s assay (typically a
  `SummarizedExperiment` or a `SingleCellExperiment`) such as
  `process(object)`, returns a new object of the same type.

- A processing function such `process(object, i)`, that acts on a
  Feautre object takes a second argument `i`, that can be a vector of
  indices or names, returns a new object of class QFeatures with its
  assay(s) `i` modified according to `process(object[[i]])`.

- The argument `i` mustn’t be missing, i.e. one shouldn’t (at least in
  general) permit to (blindly) apply some processing on all assays.

### Assays

- Assays should have unique rownames (even though this isn’t required
  for SEs). If they aren’t, only the first occurence of the name is
  kept:

&nbsp;

    hlpsms <- hlpsms[1:5000, ] ## faster

    ft1 <- readQFeatures(hlpsms, ecol = 1:10, name = "psms", fname = "Sequence")
    sum(rownames(ft1[[1]]) == "ANLPQSFQVDTSk")
    ft1 <- aggregateFeatures(ft1, "psms", fcol = "Sequence",
                             name = "peptides", fun = colSums)
    sapply(rownames(ft1), anyDuplicated)
    ft1

    ## subsetting still works
    ft2 <- subsetByFeature(ft1, "ANLPQSFQVDTSk")
    ft2

The underlying reason why this fails is due to matrix subsetting by name
when these names aren’t unique.

``` r

m <- matrix(1:10, ncol = 2)
colnames(m) <- LETTERS[1:2]
rownames(m) <- c("a", letters[1:4])
m
m["a", ]
```

And of course, this affects SEs …

``` r

se <- SummarizedExperiment(m)
assay(se["a", ])
```

… and MultiAssayExperiments.

Note that in the example above, `"ANLPQSFQVDTSk"` is present in both the
`psms` and `peptides` assays, and the

    for (k in setdiff(all_assays_names, leaf_assay_name)) { ... }

loop in `.subsetByFeature` isn’t executed at all. This will need to be
investigated. But the behaviour above can be reproduced even when that’s
not the case. See

    hlpsms$Sequence2 <- paste0(hlpsms$Sequence, "2")
    ft1 <- readQFeatures(hlpsms, ecol = 1:10, name = "psms", fname = "Sequence2")
    ...

This could be **fixed** by switching to indices:

    > (i <- which(rownames(m) == "a"))
    [1] 1 2
    > m[i, ]
      A B
    a 1 6
    a 2 7
    > se[i, ]
    class: SummarizedExperiment
    dim: 2 2
    metadata(0):
    assays(1): ''
    rownames(2): a a
    rowData names(0):
    colnames(2): A B
    colData names(0):

See issue \#91.

# Assay links

Currently, we have

- Assay links produces by `aggregateFeatures` and manually with
  `addAssayLink`.

- *One-to-one* Assay links produced by a processing function such as
  `logTransform` or with `addAssayLinkOneToOne`. These contain
  `"OneToOne"` in the `fcol` slot (issue 42).

- There will be a need for an assay link stemming from combining assays
  (see below and issue 52).

# Joining assays (closed)

To *combine* assays, we also need 1. relaxed `MatchedAssayExperiment`
constrains (see \#46) 2. assay links with multiple parent assays (see
\#52)

`combine,MSnSet,MSnSet` does two things, i.e. `rbind` and `cbind`. Here,
we nedd (at least in a first instance) and have
`cbind,SummarizedExperiment`.

- do we need some constrains requiering identical rownames?
  `cbind,SummarizedExperiment` uses the mcols to check whether rows
  match.
- should unique rows in one assay get NAs in the other one? yes!

We need a **join**-type of function, that adds NAs at the assay level.
To do this, we need to have a union of features before rbinding the
assays.

As for rowData, we want to

- keep the mcols that match exactly between assays (ex: PeptideSequence,
  ProteinAccession, …)
- remove mcols that differ between assays (ex: PEP, qvalues, charge,
  rtime, …)

The row data will be accessible through links between assays anyway.

Naming:

    joinAssays(QFeatures, c("pep_batch1", "pep_batch2", "pep_batch3"), name = "peptides")
    joinAssays(QFeatures, c(1, 2, 3), name = "peptides")

Algorithm: 1. Find which mcols to keep 2. Extend with rownames and NAs
(depending on type of join) 3. Order assays 4. cbind assays (see
`cbind,SummarizedExperiment`)

Do we want a public *join* for SummarizedExperiments? Discuss with SE
maintainers.

Note: if we were to have assay from multiple fractions to be *rbind*ed,
we could consider a `rbindAssays`, `mergeFractions`, `bindFractions`, …

# Replacing vs adding assays

Issues
[https://github.com/rformassspectrometry/QFeatures/issues/193](https://rformassspectrometry.github.io/QFeatures/193)
and
[https://github.com/rformassspectrometry/QFeatures/issues/186](https://rformassspectrometry.github.io/QFeatures/186).

Currently, assays are replaced with - filterNA() - filterFeatures() (and
possibly others)

Sometimes, we want to add, rather than replace, for example if we want
to test/assess the effect of different filters. This could be defined by
the `names` argument. If missing (default), the assays are replaced. If
present and of same length than `i`, new assays are added.

- When it comes to data processing, we could also have a *subset*
  argument, that would implicitly only process a subset of rows so as to
  avoid to explicitly store the subset/intermediate assay.

- A more radical change would be for
  [`filterFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-filtering.md)
  to add a rowData logical that defines the rows to be filtered.

There are multiple ideas/discussion replated to QFeatures becoming very
large (and slow). Rather than adding more assays, we could: - use
logical for subsetting; - use multiple assays within a
SingleCellExperiment (or SE), when the dimensions remain identical (for
exmple logTransform()); - have a unique database to handle and manage
all data (assays and rowData).

But we agree that the interface, for the user, should remain simple,
i.e. different assays. For now, keep the same philosophy and create new
assays for all operations, and start a reflexion for more in-depth
refactoring.

See also [HDF5 backend
issue](https://github.com/rformassspectrometry/QFeatures/issues/171).

# Devel roadmap

## scp

- scp data [import](https://github.com/UCLouvain-CBIO/scp/issues/48)
  documentation. Current approach doesn’t work well for LFQ. Add
  documentation, and align scpGUIImport.
- [Compute mertics per
  cell](https://github.com/UCLouvain-CBIO/scp/issues/44) is clumsy/ad
  hoc at the moment. We need a general function to iterates over
  cells/features and store automatically in col/rowData. What is
  computed defined by a user-defined function. For example
  `*OverSample()` and `*OverFeatures()` and a more general
  `computeMetric()`. These would added to `QFeatures` with specialised
  metrics implemented in `scp`.
- Merge scplainer into scp (issue
  [45](https://github.com/UCLouvain-CBIO/scp/pull/45)): needs unit tests
- [Missing data
  vignette/report](https://github.com/UCLouvain-CBIO/scp/issues/40):
  unit tests coming soon
- [readScpFromDIANN()](https://github.com/UCLouvain-CBIO/scp/issues/42)
  needs unit test and some refactoring
- [`divideByReference()`](https://github.com/UCLouvain-CBIO/scp/issues/34)
  is a batch correction method, and should be moved to QFeatures, and
  made accessible as part of a more general `batchCorrect()` interface,
  that can also be used for combat and limma.

## scpdata

- Keep [adding new
  datasets](https://github.com/UCLouvain-CBIO/scpdata/issues/15).

## QFeatures

- [All QFeatures functions should add
  assays](https://github.com/rformassspectrometry/QFeatures/issues/193) -
  see ‘Replacing vs adding assays’ above.
- [`filterNA()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  on columns (see
  [https://github.com/rformassspectrometry/QFeatures/issues/173](https://rformassspectrometry.github.io/QFeatures/173)\].
  Or should we use `computeMetricOverCells()` with metric `nNA`?
- Small bugs: [nNA() on empty
  assays](https://github.com/rformassspectrometry/QFeatures/issues/174)
  and [longFormat(se,
  rowvars)](https://github.com/rformassspectrometry/QFeatures/issues/196)
  (eventhough the latter should probbaly use tidy versions).
- Interest in
  [pseudo-bulking](https://github.com/rformassspectrometry/QFeatures/issues/188),
  which should work out-of-the-box, since the colData is unique for all
  assaya (Chris to address during post-doc).

# Tabular input (issue 199)

1.  Single-set case, multiplexed: requires `colAnnotation` only. Also LF
    with a re-ordered peptide/protein-level table (runs are missing in
    this case).

&nbsp;

    |------+------------+-----------|
    | cols | Quant 1..N | more cols |
    |      |            |           |
    |      |            |           |
    |      |            |           |
    |------+------------+-----------|

    readQFeatures(hlpsms, quantCols = 1:10)
    readQFeatures(hlpsms, colAnnotation = colann)

    ## also possible, but redundant
    readQFeatures(hlpsms, colAnnotation = colann, quantCols = 1:10)

2.  Multi-set case, multiplexed: requires `colAnnotation` and `runCol`.

&nbsp;

    |-----+------+------------+-----------|
    | Run | cols | Quant 1..N | more cols |
    |   1 |      |            |           |
    |   1 |      |            |           |
    |-----+------+------------+-----------|
    |   2 |      |            |           |
    |-----+------+------------+-----------|

    readQFeatures(hlpsms, quantCols = 1:10, runCol = "file")
    readQFeatures(hlpsms, colAnnotation = colann, runCol = "file")

3.  Multi-set case, LF: requires `colData` and `runCol` with a optional
    `multiplexing` (for plexDIA).

&nbsp;

    |-----+------+---------+-----------+-----------|
    | Run | cols | Quant 1 | more cols | multiplex |
    |   1 |      |         |           |           |
    |   1 |      |         |           |           |
    |-----+------+---------+-----------+-----------|
    |   2 |      |         |           |           |
    |-----+------+---------+-----------+-----------|

4.  Special case DIANN. A specialised function that parses the table to
    case 2.

Users can either use the arguments above or a `colAnnotation` data.frame
(that will become the `colData`).

## DIANN data

- DIANN data is in long format
- readQFeaturesFromDIANN() would need to transform from long to wide,
  then use readQFeatures()
- multiplexing: NULL or “none”

&nbsp;

    dfr |>
      diannWider() |>
      readQFeatures()

    readQFeaturesFromDIANN <- funtion(dfr, multiplexing = NULL, ...) {
        if (!is.null(multiplexing))
            x <- .diannWider(multiplexing)
        readQFeatures(x, ...)
    }

# Possible optimsations

Simple ones

- In
  [`longForm()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-longForm.md)
  have a `na.rm` argument and propagate it to
  [`reshape2::melt()`](https://rdrr.io/pkg/reshape2/man/melt.html) for
  each assay. This would require to loop over assay, extract them with
  `getWithColData()`, and use our own `longFormSE()` supporting `na.rm`.
- Parallelise
  [`aggregateFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md),
  [`logTransform()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md),
  [`scaleTransform()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md),
  [`impute()`](https://rformassspectrometry.github.io/QFeatures/reference/impute.md),
  …
