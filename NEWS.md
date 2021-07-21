# QFeatures 1.3.0

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
