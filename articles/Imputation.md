# Imputing quantitative proteomics data

Abstract

This vignette describes the multiple ways to perform imputation using
the
[`impute()`](https://rformassspectrometry.github.io/QFeatures/reference/impute.md)
methods and the underlying support function from the *MsCoreUtils*
package. This vignette is distributed under a CC BY-SA license.

## Introduction

This vignette provides a technical description of the imputation
functionality available in the *R for Mass Spectrometry* packages, in
particular
*[MsCoreUtils](https://bioconductor.org/packages/3.23/MsCoreUtils)* for
the implementation and
*[QFeatures](https://bioconductor.org/packages/3.23/QFeatures)* for the
high level application. These packages depend on other ones for the
specific imputation implementation approaches.

``` r

library(MsCoreUtils)
library(QFeatures)
```

This vignette focuses on the technical aspects of imputation, without
delving in the scientific motivations too much - see (Webb-Robertson et
al. 2015; Lazar et al. 2016; Bramer et al. 2021) for the necessary
backgroud. We will simply introduce important concepts when needed and
refer to some relevant papers for further reading.

An important concept, described among others in (Lazar et al. 2016), is
data that can be missing *at random* (MAR) or missing *not at random*
(MNAR). A MNAR feature is assumed to be missing in the data because is
was effectively absent or below the limit of detection in the biological
sample. MAR features, on the other hand, have not been detected or
identified due to technological limitations such as poor ionisation,
competition among precursors (in data dependent acquisition), or absence
of identification or mis-identification. Given the different underlying
causes of the missingness, they should be imputed using different
approaches. Typically, MNAR features can imputed using left-censored
method, that will impute using a *small* value, reflective of the
absence of the feature, while MAR features should be imputed using *hot
deck* approaches, i.e. methods that impute using *similar* or *matching*
values.

We would like to caution users on the risks of imputation, in particular
when a high proportion of values are missing. Given the different types
of missingness, wrongly imputing values can substantially distort
downstream analyses and their validity. In such situations, it might be
safer to avoid imputation altogether, and maintain missing values. Is
can also be helpful to filter out features with *too many* missing
values - the
[`filterNA()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
function can be used to such effect.

The imputation methods available in the
*[MsCoreUtils](https://bioconductor.org/packages/3.23/MsCoreUtils)*
package can be listed programmatically with the
[`imputeMethods()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
function and are documented in the
[`?impute_matrix`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
[documentation
page](https://rformassspectrometry.github.io/MsCoreUtils/reference/imputation.html).

``` r

imputeMethods()
```

    ##  [1] "bpca"    "knn"     "QRILC"   "MLE"     "MLE2"    "MinDet"  "MinProb"
    ##  [8] "min"     "zero"    "mixed"   "nbavg"   "with"    "RF"      "none"

Note that 0s are technically impossible to be recorded by a mass
spectrometer, and should never be observed in a dataset. If present,
these are the result of a prior zero-imputation by the pre-processing
software that erroneously suggest that the feature was of the MNAR type
and effectively absent in the sample. We advise to start your processing
by replacing these misleading 0 by a missing value. This could be
achieved with the
[`zeroIsNA()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
method if your data is formated as a `SummarizedExperiment` object
(Morgan et al. 2019).

### Example data

``` r

m <- matrix(1:50, nrow = 10)
diag(m) <- NA
m[which(is.na(m)) + 5] <- NA
dimnames(m) <- list(paste0("F", 1:10), paste0("S", 1:5))
randna <- rep(c(TRUE, FALSE), each = 5)
se <- SummarizedExperiment(assays = m,
                           rowData = data.frame(randna))
se
```

    ## class: SummarizedExperiment 
    ## dim: 10 5 
    ## metadata(0):
    ## assays(1): ''
    ## rownames(10): F1 F2 ... F9 F10
    ## rowData names(1): randna
    ## colnames(5): S1 S2 S3 S4 S5
    ## colData names(0):

We are going to use the small `SummarizedExperiment` to illustrate the
different imputation approaches and their parametrisation. It is
composed of 10 features and 5 samples, and contains 5 missing values
aligned diagonally along the top and bottom parts of the data matrix.

``` r

assay(se)
```

    ##     S1 S2 S3 S4 S5
    ## F1  NA 11 21 31 41
    ## F2   2 NA 22 32 42
    ## F3   3 13 NA 33 43
    ## F4   4 14 24 NA 44
    ## F5   5 15 25 35 NA
    ## F6  NA 16 26 36 46
    ## F7   7 NA 27 37 47
    ## F8   8 18 NA 38 48
    ## F9   9 19 29 NA 49
    ## F10 10 20 30 40 NA

In the following sections, we will use the [`impute()`
method](https://rformassspectrometry.github.io/QFeatures/reference/impute.html)
(see
[`?QFeatures::impute`](https://rformassspectrometry.github.io/QFeatures/reference/impute.md))
and apply it to the `SummarizedExperiment` instance `se` above. The
method is also applicable to `QFeatures` objects. The individual
[`impute_matrix()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
and other `impute_*` from the `MsCoreUtils` package can be applied
directly on `matrix` objects.

## Simple imputation

We refer to *simple* (or *single*) imputation when a single imputation
method is used across the whole dataset. For example, if we want to
replace all missing values by 0, we can use the
[`impute()`](https://rformassspectrometry.github.io/QFeatures/reference/impute.md)
method as shown below.

``` r

impute(se, method = "zero") |> assay()
```

    ##     S1 S2 S3 S4 S5
    ## F1   0 11 21 31 41
    ## F2   2  0 22 32 42
    ## F3   3 13  0 33 43
    ## F4   4 14 24  0 44
    ## F5   5 15 25 35  0
    ## F6   0 16 26 36 46
    ## F7   7  0 27 37 47
    ## F8   8 18  0 38 48
    ## F9   9 19 29  0 49
    ## F10 10 20 30 40  0

Setting the `method` to `"zero"` will apply the
[`MsCoreUtils::impute_zero()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
function on the object’s assay.

### Passing paramters to the imputation function

Or, if we want to impute all missing values with a specific value such
as 0.5, we can use the `"with"` method to apply the
[`MsCoreUtils::impute_with()`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
function. This function requires an additional argument, `val`, that
defines the specific value that should be used to replace missing
values.

``` r

impute(se, method = "with", val = 0.5) |> assay()
```

    ##       S1   S2   S3   S4   S5
    ## F1   0.5 11.0 21.0 31.0 41.0
    ## F2   2.0  0.5 22.0 32.0 42.0
    ## F3   3.0 13.0  0.5 33.0 43.0
    ## F4   4.0 14.0 24.0  0.5 44.0
    ## F5   5.0 15.0 25.0 35.0  0.5
    ## F6   0.5 16.0 26.0 36.0 46.0
    ## F7   7.0  0.5 27.0 37.0 47.0
    ## F8   8.0 18.0  0.5 38.0 48.0
    ## F9   9.0 19.0 29.0  0.5 49.0
    ## F10 10.0 20.0 30.0 40.0  0.5

Each of the underlying function’s details are documented in the
[`?impute_zero`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html),
[`?impute_with`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html), …
manual pages, that all lead to the main
[`?impute_matrix`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
documentation.

## The MARGIN argument

In the two simple examples above, there is no sense of direction when
imputing, as every missing value is replaced by a single, pre-defined
value. In many cases however, this is not the case. To illustrate this,
let’s use the `"MinDet"` method (see
[`?impute_MinDet`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)).
*MinDet* performs the imputation of left-censored missing data using a
deterministic minimal value approach. Considering a expression data with
*n* samples and *p* features, for each *sample*, the missing entries are
replaced with a minimal value observed in that sample. The minimal value
observed is estimated as being the q-th quantile (default `q = 0.01`) of
the observed values in that sample.

Below, we are going to set `q = 0` to impute with the minimal value
within each sample.

``` r

impute(se, method = "MinDet", q = 0) |> assay()
```

    ## Imputing along margin 2 (samples/columns).

    ##     S1 S2 S3 S4 S5
    ## F1   2 11 21 31 41
    ## F2   2 11 22 32 42
    ## F3   3 13 21 33 43
    ## F4   4 14 24 31 44
    ## F5   5 15 25 35 41
    ## F6   2 16 26 36 46
    ## F7   7 11 27 37 47
    ## F8   8 18 21 38 48
    ## F9   9 19 29 31 49
    ## F10 10 20 30 40 41

As can be seen, the missing values in sample S1, namely F1 and F6, have
been imputed by the smallest observed value in S1, namely 2. And
similarly for the four other samples.

In the definition above, it is explicitly stated that the imputation is
done for each sample, i.e. along the columns of the quantitative matrix,
also called the second margin. We can repeat the same imputation by
explicitly setting `MARGIN = 2`.

``` r

impute(se, method = "MinDet", q = 0, MARGIN = 2) |> assay()
```

    ## Imputing along margin 2 (samples/columns).

    ##     S1 S2 S3 S4 S5
    ## F1   2 11 21 31 41
    ## F2   2 11 22 32 42
    ## F3   3 13 21 33 43
    ## F4   4 14 24 31 44
    ## F5   5 15 25 35 41
    ## F6   2 16 26 36 46
    ## F7   7 11 27 37 47
    ## F8   8 18 21 38 48
    ## F9   9 19 29 31 49
    ## F10 10 20 30 40 41

And indeed, the default margin for the `"MinDet"` method is 2:

``` r

getImputeMargin("impute_MinDet")
```

    ## [1] 2

The imputation margin is not always 2. The *nearest neighbour*
imputation method chooses a certain number of similar features. By
similar features, we explicitly refer to other rows, i.e. the first
margin:

``` r

getImputeMargin("impute_knn")
```

    ## [1] 1

It is possible to change the margin from its default value. Below, we
now use `"MinDet"` and choose the smallest value within each
feature/row.

``` r

impute(se, method = "MinDet", q = 0, MARGIN = 1) |> assay()
```

    ## Imputing along margin 1 (features/rows).

    ##     S1 S2 S3 S4 S5
    ## F1  11 11 21 31 41
    ## F2   2  2 22 32 42
    ## F3   3 13  3 33 43
    ## F4   4 14 24  4 44
    ## F5   5 15 25 35  5
    ## F6  16 16 26 36 46
    ## F7   7  7 27 37 47
    ## F8   8 18  8 38 48
    ## F9   9 19 29  9 49
    ## F10 10 20 30 40 10

Now, we see that the missing F1 value in S1 has been imputed by the
smallest observed value along the first row, namely 11.

We can extract all default margin values for all `MsCoreUtils::impute_*`
functions as show below.

``` r

getImputeMargin()
```

    ## $impute_bpca
    ## [1] 1
    ## 
    ## $impute_fun
    ## [1] 1
    ## 
    ## $impute_knn
    ## [1] 1
    ## 
    ## $impute_matrix
    ## [1] NA
    ## 
    ## $impute_min
    ## [1] NA
    ## 
    ## $impute_MinDet
    ## [1] 2
    ## 
    ## $impute_MinProb
    ## [1] 2
    ## 
    ## $impute_mixed
    ## c(1L, 1L)
    ## 
    ## $impute_mle
    ## [1] 2
    ## 
    ## $impute_neighbour_average
    ## [1] 1
    ## 
    ## $impute_QRILC
    ## [1] 2
    ## 
    ## $impute_RF
    ## [1] 2
    ## 
    ## $impute_with
    ## [1] NA
    ## 
    ## $impute_zero
    ## [1] NA

A missing margin means that, as for `"with"` or `"zero"` above, there is
not margin along which the imputation is performed. Mixed imputation is
a special case that has two margins, which we will describe in the next
section.

The relevance of the imputation margin can also depend on downstream
analyses. In (Vanderaa and Gatto 2023), the authors illustrate that
imputation along the first margin increases the correlation between
features, while imputation along the second margin increases the
correlation between samples. These artificially improved correlations
can then in turn impact any analyses that rely on the identification of
sample or protein clusters.

## Mixed imputation

As we have seen above, different underlying processes can lead to
different types of missing values, namely MAR and MNAR. One view of this
is to define these processes at the feature level [^1]. In such cases,
one might want to impute different sets of features in a *mixed* way:
MAR features with a MAR-appropriate method, and MNAR features with a
MNAR-appropriate method. This is possible is the `"mixed"` method.

To be able to apply mixed imputation, we need to define features that
are MAR, and features that are MNAR. This is done using a logical vector
whose length is equal to the number of features.

``` r

rowData(se)$randna
```

    ##  [1]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE

The `TRUE` values define the MAR features (F1 to F5 in our case) and
`FALSE` defines MNAR features (F6 to F10).

To use mixed imputation, we need to specify the MAR and MNAR features,
two imputation methods, one for MAR features, and another one for MNAR
features.

``` r

impute(se, method = "mixed",
       randna = rowData(se)$randna,
       mar = "MinDet",
       mnar = "zero") |>
    assay()
```

    ## Imputing along margin 1 (features/rows).

    ##       S1   S2   S3   S4   S5
    ## F1  11.3 11.0 21.0 31.0 41.0
    ## F2   2.0  2.6 22.0 32.0 42.0
    ## F3   3.0 13.0  3.3 33.0 43.0
    ## F4   4.0 14.0 24.0  4.3 44.0
    ## F5   5.0 15.0 25.0 35.0  5.3
    ## F6   0.0 16.0 26.0 36.0 46.0
    ## F7   7.0  0.0 27.0 37.0 47.0
    ## F8   8.0 18.0  0.0 38.0 48.0
    ## F9   9.0 19.0 29.0  0.0 49.0
    ## F10 10.0 20.0 30.0 40.0  0.0

We can see that the bottom-half of the matrix corresponding to MNAR
features have been imputed by zero, while the other have imputed by
`"MinDet"`. We also see that the margin used for `"MinDet"` was 1 (along
the rows). Indeed, the default margins are 1 for both MAR and MNAR
features: the first one is for MAR, and the second one for MNAR.

``` r

getImputeMargin("impute_mixed")
```

    ## c(1L, 1L)

### Different margins

It is of course possible to change the margins when performing mixed
imputation:

``` r

impute(se, method = "mixed",
       randna = rowData(se)$randna,
       mar = "MinDet",
       mnar = "zero",
       MARGIN = c(2, NA)) |>
    assay()
```

    ## Imputing along margin 2 (samples/columns).

    ##        S1    S2    S3    S4    S5
    ## F1   2.03 11.00 21.00 31.00 41.00
    ## F2   2.00 11.06 22.00 32.00 42.00
    ## F3   3.00 13.00 21.03 33.00 43.00
    ## F4   4.00 14.00 24.00 31.03 44.00
    ## F5   5.00 15.00 25.00 35.00 41.03
    ## F6   0.00 16.00 26.00 36.00 46.00
    ## F7   7.00  0.00 27.00 37.00 47.00
    ## F8   8.00 18.00  0.00 38.00 48.00
    ## F9   9.00 19.00 29.00  0.00 49.00
    ## F10 10.00 20.00 30.00 40.00  0.00

We set `NA` for zero-imputation (it could also have been 1, as it is
irrelevant anyway) and 2 for MinDet-imputation. And we can confirm that
this time, the the MAR features have been imputed using the smallest
values have been choosen for each sample/column.

### Passing paramters to the imputation functions

It is possible to pass arguments to the respective MAR and MNAR function
using the `marArgs` and `mnarArg` arguments as named lists. Below, we
are going to use *MinDet* in both cases, with different parameters.

``` r

impute(se,
       method = "mixed",
       randna = rowData(se)$randna,
       mar = "MinDet",
       mnar = "MinDet",
       marArgs = list(q = 0),
       mnarArgs = list(q = 1),
       MARGIN = c(1, 1)) |>
    assay()
```

    ## Imputing along margin 1 (features/rows).
    ## Imputing along margin 1 (features/rows).

    ##     S1 S2 S3 S4 S5
    ## F1  11 11 21 31 41
    ## F2   2  2 22 32 42
    ## F3   3 13  3 33 43
    ## F4   4 14 24  4 44
    ## F5   5 15 25 35  5
    ## F6  46 16 26 36 46
    ## F7   7 47 27 37 47
    ## F8   8 18 48 38 48
    ## F9   9 19 29 49 49
    ## F10 10 20 30 40 40

In both cases, we impute along the rows. For the MAR features (top half
of the matrix), we impute using the minimal value of that row (using
`q = 0`), while for the MNAR feature (bottom half of the matrix), we
impute using the maximal value of that row (using `q = 1`). As
anticipated, the value of F1 in S1 gets 11, and F5 in S1 gets 46.

### Using the whole matrix to compute imputated values

When doing mixed imputation, the respective MAR and MNAR sub-matrices
are split and imputed separately. It is also possible the use the whole
data matrix to compute the MAR and MNAR imputated values. This is
controlled by the `split` argument that, by default, is set to `TRUE`.

Below, we are going to repeat a mixed imputation, imputing the MAR
values (the top half of the matrix) with the highest value of the
*whole* columns using `MARGIN = 2` and `split = TRUE`. The NMAR values
(the bottom half of the matrix) are impute using the smallest value
along the rows using `MARGIN = 1`, and are hence not impacted by the
`split` value.

``` r

impute(se,
       method = "mixed",
       randna = rowData(se)$randna,
       mar = "MinDet",
       mnar = "MinDet",
       marArgs = list(q = 1),
       mnarArgs = list(q = 0),
       MARGIN = c(2, 1),
       split = FALSE) |>
    assay()
```

    ## Imputing along margin 2 (samples/columns).

    ## Imputing along margin 1 (features/rows).

    ##     S1 S2 S3 S4 S5
    ## F1  10 11 21 31 41
    ## F2   2 20 22 32 42
    ## F3   3 13 30 33 43
    ## F4   4 14 24 40 44
    ## F5   5 15 25 35 49
    ## F6  16 16 26 36 46
    ## F7   7  7 27 37 47
    ## F8   8 18  8 38 48
    ## F9   9 19 29  9 49
    ## F10 10 20 30 40 10

We see that the value of F1 in S1 gets 10, the highest value from F10 in
S1. If we keep the default `split = TRUE`, it would have gotten 5 from
F5, the highest value among the MAR values. The MNAR imputation isn’t
affected by the split and get the smallest values in each row.

## Session information

    ## R version 4.6.0 (2026-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] QFeatures_1.23.1            MultiAssayExperiment_1.38.0
    ##  [3] SummarizedExperiment_1.42.0 Biobase_2.72.0             
    ##  [5] GenomicRanges_1.64.0        Seqinfo_1.2.0              
    ##  [7] IRanges_2.46.0              S4Vectors_0.50.0           
    ##  [9] BiocGenerics_0.58.0         generics_0.1.4             
    ## [11] MatrixGenerics_1.24.0       matrixStats_1.5.0          
    ## [13] MsCoreUtils_1.25.3          BiocStyle_2.40.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyr_1.3.2             sass_0.4.10             SparseArray_1.12.2     
    ##  [4] stringi_1.8.7           lattice_0.22-9          magrittr_2.0.5         
    ##  [7] digest_0.6.39           evaluate_1.0.5          grid_4.6.0             
    ## [10] bookdown_0.46           fastmap_1.2.0           plyr_1.8.9             
    ## [13] jsonlite_2.0.0          Matrix_1.7-5            ProtGenerics_1.44.0    
    ## [16] BiocManager_1.30.27     purrr_1.2.2             lazyeval_0.2.3         
    ## [19] textshaping_1.0.5       jquerylib_0.1.4         abind_1.4-8            
    ## [22] cli_3.6.6               rlang_1.2.0             XVector_0.52.0         
    ## [25] cachem_1.1.0            DelayedArray_0.38.1     yaml_2.3.12            
    ## [28] otel_0.2.0              S4Arrays_1.12.0         tools_4.6.0            
    ## [31] reshape2_1.4.5          dplyr_1.2.1             vctrs_0.7.3            
    ## [34] R6_2.6.1                lifecycle_1.0.5         stringr_1.6.0          
    ## [37] fs_2.1.0                htmlwidgets_1.6.4       clue_0.3-68            
    ## [40] MASS_7.3-65             ragg_1.5.2              cluster_2.1.8.2        
    ## [43] pkgconfig_2.0.3         desc_1.4.3              pillar_1.11.1          
    ## [46] pkgdown_2.2.0.9000      bslib_0.10.0            Rcpp_1.1.1-1.1         
    ## [49] glue_1.8.1              systemfonts_1.3.2       tidyselect_1.2.1       
    ## [52] tibble_3.3.1            xfun_0.57               knitr_1.51             
    ## [55] AnnotationFilter_1.36.0 igraph_2.3.0            htmltools_0.5.9        
    ## [58] rmarkdown_2.31          compiler_4.6.0

## References

Bramer, Lisa M, Jan Irvahn, Paul D Piehowski, Karin D Rodland, and
Bobbie-Jo M Webb-Robertson. 2021. “A Review of Imputation Strategies for
Isobaric Labeling-Based Shotgun Proteomics.” *J. Proteome Res.* 20 (1):
1–13.

Lazar, Cosmin, Laurent Gatto, Myriam Ferro, Christophe Bruley, and
Thomas Burger. 2016. “Accounting for the Multiple Natures of Missing
Values in Label-Free Quantitative Proteomics Data Sets to Compare
Imputation Strategies.” *J. Proteome Res.* 15 (4): 1116–25.

Morgan, Martin, Valerie Obenchain, Jim Hester, and Hervé Pagès. 2019.
*SummarizedExperiment: SummarizedExperiment Container*.

Vanderaa, Christophe, and Laurent Gatto. 2023. “Revisiting the Thorny
Issue of Missing Values in Single-Cell Proteomics.” *arXiv
\[q-Bio.QM\]*, April.

Webb-Robertson, Bobbie-Jo M, Holli K Wiberg, Melissa M Matzke, et al.
2015. “Review, Evaluation, and Discussion of the Challenges of Missing
Value Imputation for Mass Spectrometry-Based Label-Free Global
Proteomics.” *J. Proteome Res.* 14 (5): 1993–2001.

[^1]: Nowadays, I believe that this feature-level representation of MAR
    or MNAR is not entirely correct. It can be useful in some simple
    cases though.
