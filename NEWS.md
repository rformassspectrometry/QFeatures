# QFeatures 1.11

## QFeatures 1.11.2

- Update message to fix test upon recent changes in MAE.

## QFeatures 1.11.1

- Update `nNA()` and `filterNA()` and man pages to clarify percentages
  and proportions (see #189).

## QFeatures 1.11.0

- New Bioconductor 3.18 (stable) release

# QFeatures 1.10

## QFeatures 1.10.0

- New Bioconductor 3.17 (stable) release

# QFeatures 1.9

## QFeatures 1.9.4

- New `dropEmptyAssays()` function (see issue #184).

## QFeatures 1.9.3

- Minor rephrasing in vignette and README.

## QFeatures 1.9.2

- feat: filterFeatures() now allows to select assays to filter (i
  argument)
- feat: aggregateFeatures() can now take multiple assays
- feat: impute() can now take multiple assays
- feat: processing functions (normalize, scaleTransform, logTransform,
  sweep) can now take multiple assays
- refactor: avoid validObject() when possible
- Use `|>` rather than `%>%`.

## QFeatures 1.9.1

- fix: solved bug in `selectRowData()`

## QFeatures 1.9.0

- New Bioconductor 3.17 (devel) release

# QFeatures 1.8

## QFeatures 1.8.0

- New Bioconductor 3.16 (stable) release

# QFeatures 1.7

## QFeatures 1.7.3

- fix: fixed filterFeatures when selection contains environment variables

## QFeatures 1.7.2

- feat: added `c` methods to combine QFeatures objects.
- feat: added `nrows` and `ncols` methods. Also added use.names argument
  (cf ?BiocGenerics::dims)
- docs: improved docs for `filterFeatures()`
- tests: improved unit tests for `filterFeatures()`
- feat: added a keep argument in `filterFeatures()` to control whether
  to keep or remove features for assays that do not contain the filter
  variable. Also added message printing for a better overview of which
  variable were found.
- fix: fixed `addAssay()` to solve issue #104.
- refactor: refactored `addAssay()` and dramatically improved the usage
  of computational resources.
- feat: `colData` is automatically transferred from the assay to the
  QFeatures object.
- feat: implemented `removeAssay()` and `replaceAssay()`. Together with
  `addAssay()`, these functions are used to implement the replacement
  method `[[<-` required to solve issue #57.
- Add CC-BY-SA license for vignettes.

## QFeatures 1.7.1

- refactor: imputation now adds a new assay instead of replacing values.

## QFeatures 1.7.0

- New Bioc devel version.

# QFeatures 1.5

## QFeatures 1.5.3

- feat: aggregation by adjacency matrix
- New `adjacencyMatrix,SummarizedExperiment` and
  `adjacencyMatrix,QFeatures` methods using
  `ProtGenerics::adjacencyMatrix` (available in version >= 1.27.1).

## QFeatures 1.5.2

- fix: implemented an `updateObject()` method for `QFeatures` objects.

## QFeatures 1.5.1

- Document the use of peptide/protein adjacency matrices in
  `aggregateFeatures()` and new `adjacencyMatrix()` accessor.

## QFeatures 1.5.0

- New devel version (Bioc 3.15)

# QFeatures 1.3

## QFeatures 1.3.6

- New `feat3` example data to demonstrate and test more complex AssayLinks
  structure.
- Improved the `plot,QFeautres` function to avoid cluttering of nodes.
- Adapted the visualization vignette using `feat3`.

## QFeatures 1.3.5

- Add plot,QFeatures and visualisation vignette.

## QFeatures 1.3.4

- Fixed bug that produced invalid AssayLinks when using filterNA.

## QFeatures 1.3.3

- Improved validity checks on `AssayLinks`
- Fixed the subsetting of `AssayLinks` to ensure consistent data

## QFeatures 1.3.2

- Add logo to package
- Fix class coercion error (see #b9ce7f1e9)

## QFeatures 1.3.1

- Added `rbindRowData`: a function to select variables in the `rowData`
  and bind it in a single `DataFrame`
- Added `rowData<-`: this new method replaces `replaceRowDataCols` to
  offer a more standardize functionality.
- Added a new section in the `QFeatures` vignette to expand on how to
  manipulate the metadata within a `QFeatures` object

## QFeatures 1.3.0

- New devel version (Bioc 3.14)

# QFeatures 1.2.0

## QFeatures 1.2.0

- New release version (Bioc 3.13)

# QFeatures 1.1.0

## QFeatures 1.1.4

- Added `replaceRowDataCols` and `removeRowDataCols`, two functions
  to streamline manipulation of rowData within a `QFeature` object.

## QFeatures 1.1.3

- Added `countUniqueFeatures`, a function to count the number of
  unique features per sample.

## QFeatures 1.1.2

- Manually install preprocessCore (see
  https://github.com/Bioconductor/bioconductor_docker/issues/22 for
  details) to use quantile normalisation in vignette and tests.

- Update vignette to show `normalize()` and `logTransform()` directly
  on a `QFeatures` object and reference the QFeaturesWorkshop2020
  workshop and WSBIM2122 chap 8.

## QFeatures 1.1.1

- Fix typo in vignette
- Improve first paragraph in intro vignette.

## QFeatures 1.1.0

- New Bioconductor devel version

# QFeatures 0.99

## QFeatures 0.99.4

- Fix: improved `nNA` with new implementation and additional unit
  tests <2020-10-23 Fri>

## QFeatures 0.99.3

- New feature: the `longFormat` function returns a long `DataFrame`
  with quantitative data along with metadata (see #116)
  <2020-10-8 Thu>
- New feature: the `rowData` method returns a list containing the
  `rowData` for all assays (see #86)
  <2020-09-16 Wed>
- Keep colnames when reading a single column assay (see #108)
  <2020-09-09 Wed>

## QFeatures 0.99.2

- Added `infIsNA`.
- Add Christophe Vanderaa as an author.

## QFeatures 0.99.1

- Address comments Bioconductor review (see [submission
  issue](https://github.com/Bioconductor/Contributions/issues/1556)
  for details.

## QFeatures 0.99.0

- Bioconductor submission
