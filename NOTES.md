# One or two levels in Features (closed)

One or two assay levels could be considered in Features:

- one level: each SE contains only a single assay, and when an SE is
  processed (log-transformed, normalised, ...) in a way that its
  dimensions stay the same, a new SE is created and added to the
  Features object.
  
- two level: SEs can contain multiple assays, and when an SE is
  processed (log-transformed, normalised, ...) in a way that its
  dimensions stay the same, a new assay is added to that SE.
  
This
[question](https://stat.ethz.ch/pipermail/bioc-devel/2020-January/016096.html)
on the bioc-devel list ask for advice on SE processing, and whether a
new SE or new assay in the original SE should be preferred. While the
letter is arguably more elegant, and is also used in
SingleAssayExperiment pipelines, it doesn't seem to be the case when
using SummarizedExperiments.

As for features (or [MultiAssayExperiments in
general](https://github.com/waldronlab/MultiAssayExperiment/issues/266)),
the two-level approach isn't readily available out-of-the-box, and
would require additional developments:

- Every function that operates on an SE of a Features object would
  need to allow the user to specify which assay to use (and/or by
  default use the latest one).
  
- The `show,Features` method would need to display the number/names of
  the assays in each SE to make these two levels explicit.

Despite the elegant of the two-level option, it seems that the
additional development isn't warranted at this time.

The [`updateAssay`
function](https://github.com/rformassspectrometry/Features/issues/37)
was originally intended for the two-level approach, i.e. to add an
assay to an SE. This is not considered anymore (for now, at least).

# How to add new assays

1. Through aggregation with `aggregateFeatures`.

2. Processing an SE. 

This can/could be done explicitly with `addAssay`

```
addAssay(cptac, logTransform(cptac[["peptides"]]), name = "peptides_log")
addAssay(cptac, logTransform(cptac[[1]]), name = "peptides_log")
```

or implicitly 

```
logTransform(cptac, "peptides", name = "peptides_log")
logTransform(cptac, 1, name = "peptides_log")
```

3. Joining SEs (for example multiple TMT batches) (TODO) 

```
joinAssays(Features, c("pep_batch1", "pep_batch2", "pep_batch3"), name = "peptides")
joinAssays(Features, c(1, 2, 3), name = "peptides")
```

See below.

# Features API

### Processing functions

- A processing function that acts on a Feature's assay (typically a
  `SummarizedExperiment` or a `SingleCellExperiment`) such as
  `process(object)`, returns a new object of the same type.
  
- A processing function such `process(object, i)`, that acts on a
  Feautre object takes a second argument `i`, that can be a vector of
  indices or names, returns a new object of class Features with its
  assay(s) `i` modified according to `process(object[[i]])`.
  
- The argument `i` mustn't be missing, i.e. one shouldn't (at least in
  general) permit to (blindly) apply some processing on all assays.
  
### Assays 

- Assays should have **unique rownames** (even though this isn't
  required for SEs). This is a requirement to be able to easily and
  transparently keep links between different.


# Assay links

Currently, we have

- Assay links produces by `aggregateFeatures` and manually with
  `addAssayLink`.
  
- *One-to-one* Assay links produced by a processing function such as
  `logTransform` or with `addAssayLinkOneToOne`. These contain
  `"OneToOne"` in the `fcol` slot (isseu 42).

- There will be a need for an assay link stemming from combining
  assays (see above and issue 52).

# Joining assays

To *combine* assays, we also need 
1. relaxed `MatchedAssayExperiment` constrains (see #46)
2. assay links with multiple parent assays (see #52)

`combine,MSnSet,MSnSet` does two things, i.e. `rbind` and
`cbind`. Here, we nedd (at least in a first instance) and have
`cbind,SummarizedExperiment`.

- do we need some constrains requiering identical rownames?
  `cbind,SummarizedExperiment` uses the mcols to check whether rows
  match.
- should unique rows in one assay get NAs in the other one? yes!

We need a **join**-type of function, that adds NAs at the assay
level. To do this, we need to have a union of features before rbinding
the assays.

As for rowData, we want to 

- keep the mcols that match exactly between assays (ex:
  PeptideSequence, ProteinAccession, ...)
- remove mcols that differ between assays (ex: PEP, qvalues, charge,
  rtime, ...)

The row data will be accessible through links between assays anyway.

Naming:

```
joinAssays(Features, c("pep_batch1", "pep_batch2", "pep_batch3"), name = "peptides", type = "full")
joinAssays(Features, c(1, 2, 3), name = "peptides", type = "inner")
```

Algorithm:
1. Find which mcols to keep
2. Extend with rownames and NAs (depending on type of join)
3. Order assays
4. cbind assays (see `cbind,SummarizedExperiment`)

Do we want a public *join* for SummarizedExperiments? Discuss with SE
maintainers.

Note: if we were to have assay from multiple fractions to be
*rbind*ed, we could consider a `rbindAssays`, `mergeFractions`,
`bindFractions`, ...
