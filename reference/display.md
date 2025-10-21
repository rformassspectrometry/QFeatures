# Interactive MultiAssayExperiment Explorer

A shiny app to browser and explore the assays in an
`MultiAssayExperiment` object. Each assay can be selected from the
dropdown meny in the side panel, and the quantitative data and row
metadata are displayed in the respective *Assay* and *Row data* tabs.
The *Heatmap* tab displays a heatmap of the assay. The selection of rows
in the *Row data* table is used to subset the features displayed in the
*Assay* table and the heatmap to those currectly selected. See
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
for an example.

## Usage

``` r
display(object, n = 100, ...)
```

## Arguments

- object:

  An instance inheriting from `MultiAssayExperiment`.

- n:

  A `numeric(1)` indicating the maximum number of features (rows) to
  consider before disabling row clustering and displaying feature names
  for speed purposes. Default is 100.

- ...:

  Additional parameters (other than `Rowv` and `labRow`, which are set
  internally based on the value of `n`) passed to heatmap.

## Value

Used for its side effect.

## Author

Laurent Gatto

## Examples

``` r
if (FALSE) { # \dontrun{
data(feat2)
display(feat2)
} # }
```
