# Example `QFeatures` object after processing

`feat3` is a small `QFeatures` object that contains 7 assays: `psms1`,
`psms2`, `psmsall`, `peptides`, `proteins`, `normpeptides`,
`normproteins`. The dataset contains example data that could be obtained
after running a simple processing pipeline. You can use it to get your
hands on manipulating `AssayLinks` since all 3 general cases are
present:

- One parent to one child `AssayLink`: the relationship can either be
  one row to one row (e.g. "peptides" to "normpeptides") or multiple
  rows to one row (e.g. "peptides" to "proteins").

- One parent to multiple children `AssayLink`: for instance "peptides"
  to "normpeptides" and "proteins".

- Multiple parents to one child `AssayLink`: links the rows between
  multiple assays to a single assays where some rows in different parent
  assays may point to the same row in the child assay. E.g. "psms1" and
  "psms2" to "psmsall"

## Usage

``` r
feat3
```

## Format

An object of class `QFeatures` of length 7.

## Source

`feat3` was built from `feat1`. The source code is available in
[`inst/scripts/test_data.R`](https://github.com/rformassspectrometry/QFeatures/blob/master/inst/scripts/test_data.R)

## See also

See
[`?feat1`](https://rformassspectrometry.github.io/QFeatures/reference/feat1.md)
for other example/test data sets.

## Examples

``` r

data("feat3")
plot(feat3)

```
