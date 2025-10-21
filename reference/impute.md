# Quantitative proteomics data imputation

The `impute` method performs data imputation on `QFeatures` and
`SummarizedExperiment` instance using a variety of methods.

Users should proceed with care when imputing data and take precautions
to assure that the imputation produce valid results, in particular with
naive imputations such as replacing missing values with 0.

See
[`MsCoreUtils::impute_matrix()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
for details on the different imputation methods available and
strategies.

## Usage

``` r
impute

# S4 method for class 'SummarizedExperiment'
impute(object, method, ...)

# S4 method for class 'QFeatures'
impute(object, method, ..., i, name = "imputedAssay")
```

## Format

An object of class `standardGeneric` of length 1.

## Arguments

- object:

  A `SummarizedExperiment` or `QFeatures` object with missing values to
  be imputed.

- method:

  `character(1)` defining the imputation method. See `imputeMethods()`
  for available ones. See
  [`MsCoreUtils::impute_matrix()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
  for details.

- ...:

  Additional parameters passed to the inner imputation function. See
  [`MsCoreUtils::impute_matrix()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
  for details.

- i:

  A `logical(1)` or a `character(1)` that defines which element of the
  `QFeatures` instance to impute. It cannot be missing and must be of
  length one.

- name:

  A `character(1)` naming the new assay name. Default is `imputedAssay`.

## Examples

``` r
MsCoreUtils::imputeMethods()
#>  [1] "bpca"    "knn"     "QRILC"   "MLE"     "MLE2"    "MinDet"  "MinProb"
#>  [8] "min"     "zero"    "mixed"   "nbavg"   "with"    "RF"      "none"   

data(se_na2)
## table of missing values along the rows (proteins)
table(rowData(se_na2)$nNA)
#> 
#>   0   1   2   3   4   8   9  10 
#> 301 247  91  13   2  23  10   2 
## table of missing values along the columns (samples)
colData(se_na2)$nNA
#>  [1] 34 45 56 39 47 52 49 61 41 42 55 45 51 43 57 53

## non-random missing values
notna <- which(!rowData(se_na2)$randna)
length(notna)
#> [1] 35
notna
#>  [1]   6  20  79  88 130 187 227 231 238 264 275 317 324 363 373 382 409 437 445
#> [20] 453 456 474 484 485 492 514 516 546 568 580 594 631 648 664 671

impute(se_na2, method = "min")
#> class: SummarizedExperiment 
#> dim: 689 16 
#> metadata(3): MSnbaseFiles MSnbaseProcessing MSnbaseVersion
#> assays(1): ''
#> rownames(689): AT1G09210 AT1G21750 ... AT4G11150 AT4G39080
#> rowData names(2): nNA randna
#> colnames(16): M1F1A M1F4A ... M2F8B M2F11B
#> colData names(1): nNA

if (require("imputeLCMD")) {
  impute(se_na2, method = "QRILC")
  impute(se_na2, method = "MinDet")
}
#> Loading required package: imputeLCMD
#> Loading required package: tmvtnorm
#> Loading required package: mvtnorm
#> Loading required package: Matrix
#> 
#> Attaching package: ‘Matrix’
#> The following object is masked from ‘package:S4Vectors’:
#> 
#>     expand
#> Loading required package: gmm
#> Loading required package: sandwich
#> 
#> Attaching package: ‘sandwich’
#> The following object is masked from ‘package:generics’:
#> 
#>     estfun
#> Loading required package: norm
#> This package has some major limitations
#> (for example, it does not work reliably when
#> the number of variables exceeds 30),
#> and has been superseded by the norm2 package.
#> Loading required package: pcaMethods
#> 
#> Attaching package: ‘pcaMethods’
#> The following object is masked from ‘package:stats’:
#> 
#>     loadings
#> Loading required package: impute
#> Imputing along margin 2 (samples/columns).
#> Imputing along margin 2 (samples/columns).
#> class: SummarizedExperiment 
#> dim: 689 16 
#> metadata(3): MSnbaseFiles MSnbaseProcessing MSnbaseVersion
#> assays(1): ''
#> rownames(689): AT1G09210 AT1G21750 ... AT4G11150 AT4G39080
#> rowData names(2): nNA randna
#> colnames(16): M1F1A M1F4A ... M2F8B M2F11B
#> colData names(1): nNA

if (require("norm"))
  impute(se_na2, method = "MLE")
#> Imputing along margin 2 (samples/columns).
#> Warning: NAs introduced by coercion to integer range
#> Iterations of EM: 
#> 1...2...
#> class: SummarizedExperiment 
#> dim: 689 16 
#> metadata(3): MSnbaseFiles MSnbaseProcessing MSnbaseVersion
#> assays(1): ''
#> rownames(689): AT1G09210 AT1G21750 ... AT4G11150 AT4G39080
#> rowData names(2): nNA randna
#> colnames(16): M1F1A M1F4A ... M2F8B M2F11B
#> colData names(1): nNA

impute(se_na2, method = "mixed",
       randna = rowData(se_na2)$randna,
       mar = "knn", mnar = "QRILC")
#> Imputing along margin 1 (features/rows).
#> Imputing along margin 1 (features/rows).
#> class: SummarizedExperiment 
#> dim: 689 16 
#> metadata(3): MSnbaseFiles MSnbaseProcessing MSnbaseVersion
#> assays(1): ''
#> rownames(689): AT1G09210 AT1G21750 ... AT4G11150 AT4G39080
#> rowData names(2): nNA randna
#> colnames(16): M1F1A M1F4A ... M2F8B M2F11B
#> colData names(1): nNA

## neighbour averaging
x <- se_na2[1:4, 1:6]
assay(x)[1, 1] <- NA ## min value
assay(x)[2, 3] <- NA ## average
assay(x)[3, 1:2] <- NA ## min value and average
## 4th row: no imputation
assay(x)
#>              M1F1A    M1F4A   M1F7A  M1F11A    M1F2B    M1F5B
#> AT1G09210       NA 0.275500 0.21600 0.18525 0.465667 0.199667
#> AT1G21750 0.332000 0.279667      NA 0.16600 0.451500 0.200375
#> AT1G51760       NA       NA 0.16825 0.18825 0.459750 0.214500
#> AT1G56340 0.336733       NA      NA      NA 0.487167 0.201833

assay(impute(x, "nbavg"))
#> Assuming values are ordered.
#> Imputing along margin 1 (features/rows).
#>              M1F1A    M1F4A     M1F7A  M1F11A    M1F2B    M1F5B
#> AT1G09210 0.166000 0.275500 0.2160000 0.18525 0.465667 0.199667
#> AT1G21750 0.332000 0.279667 0.2228335 0.16600 0.451500 0.200375
#> AT1G51760 0.166000 0.167125 0.1682500 0.18825 0.459750 0.214500
#> AT1G56340 0.336733       NA        NA      NA 0.487167 0.201833
```
