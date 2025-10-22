## Quantitative features for mass spectrometry data

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/QFeatures/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/QFeatures/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](https://codecov.io/github/rformassspectrometry/QFeatures/coverage.svg?branch=master)](https://codecov.io/github/rformassspectrometry/QFeatures?branch=master)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

![](reference/figures/logo.png)

### What is QFeatures?

`QFeatures` is a [Bioconductor
package](http://bioconductor.org/packages/QFeatures) that provides
infrastructure to management and process quantitative features for
high-throughput mass spectrometry-based proteomics assays. It provides a
familiar Bioconductor user experience to manage quantitative data across
different assay levels (such as precursors, peptide spectrum matches,
peptides and proteins or protein groups) in a coherent and tractable
format.

If you are familiar with the `MSnbase` package, `QFeatures` could be
summarised with:

> Evolving `MSnSet` data towards `SummarizedExperiment`, but for
> proteomics data.

### Getting started

The `QFeatures` class is used to manage and process quantitative
features for high-throughput mass spectrometry assays. See the
[QFeatures
introduction](https://rformassspectrometry.github.io/QFeatures/articles/QFeatures.html)
to get started and the [Processing quantitative proteomics data with
QFeatures](https://rformassspectrometry.github.io/QFeatures/articles/Processing.html)
vignette for a real-life application. Visualisation of quantitative mass
spectrometry data contained in a `QFeatures` object is illustrated in
the [Data
visualisation](https://rformassspectrometry.github.io/QFeatures/articles/Visualization.html)
vignette.

### License

The `QFeatures` code is provided under a permissive [Artistic 2.0
license](https://opensource.org/licenses/Artistic-2.0). The
documentation, including the manual pages and the vignettes, are
distributed under a [CC BY-SA
license](https://creativecommons.org/licenses/by-sa/4.0/).

## Contributions

Contributions are welcome, and should ideally be provided through a
Github pull request. Feel free to discuss any more non-trivial
suggestions or changes first in an issue.

### Contributors

- [Alexey Stukalov](https://github.com/alyst): [readQFeatures(): don’t
  ignore
  colData\$quantCols](https://github.com/rformassspectrometry/QFeatures/pull/234).
