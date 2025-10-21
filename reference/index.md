# Package index

## QFeatures class

The QFeatures class

- [`QFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`show(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`plot(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`` `[`( ``*`<QFeatures>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`` `[`( ``*`<QFeatures>`*`,`*`<character>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`c(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`dims(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`nrows(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`ncols(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`rowData(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`` `rowData<-`( ``*`<QFeatures>`*`,`*`<DataFrameList>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`` `rowData<-`( ``*`<QFeatures>`*`,`*`<ANY>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`rbindRowData()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`selectRowData()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`rowDataNames()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`` `names<-`( ``*`<QFeatures>`*`,`*`<character>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`addAssay()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`removeAssay()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`replaceAssay()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`` `[[<-`( ``*`<QFeatures>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`updateObject(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  [`dropEmptyAssays()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  : Quantitative MS QFeatures
- [`longForm(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-longForm.md)
  [`longForm(`*`<SummarizedExperiment>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-longForm.md)
  : Reshape into a long data format
- [`joinAssays()`](https://rformassspectrometry.github.io/QFeatures/reference/joinAssays.md)
  : Join assays in a QFeatures object
- [`show(`*`<AssayLink>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`updateObject(`*`<AssayLinks>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`updateObject(`*`<AssayLink>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`AssayLink()`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`AssayLinks()`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`assayLink()`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`assayLinks()`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`` `[`( ``*`<AssayLink>`*`,`*`<character>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`` `[`( ``*`<AssayLinks>`*`,`*`<list>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`addAssayLink()`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  [`addAssayLinkOneToOne()`](https://rformassspectrometry.github.io/QFeatures/reference/AssayLinks.md)
  : Links between Assays
- [`countUniqueFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/countUniqueFeatures.md)
  : Count Unique Features
- [`display()`](https://rformassspectrometry.github.io/QFeatures/reference/display.md)
  : Interactive MultiAssayExperiment Explorer

## Data import

Create QFeatures and SummarizedExperiment object

- [`readSummarizedExperiment()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.md)
  [`readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.md)
  : QFeatures from tabular data
- [`readQFeaturesFromDIANN()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeaturesFromDIANN.md)
  : Read DIA-NN output as a QFeatures objects

## Missing data

Filtering and handling of missing values

- [`zeroIsNA(`*`<SummarizedExperiment>`*`,`*`<missing>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`zeroIsNA(`*`<QFeatures>`*`,`*`<integer>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`zeroIsNA(`*`<QFeatures>`*`,`*`<numeric>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`zeroIsNA(`*`<QFeatures>`*`,`*`<character>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`infIsNA(`*`<SummarizedExperiment>`*`,`*`<missing>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`infIsNA(`*`<QFeatures>`*`,`*`<integer>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`infIsNA(`*`<QFeatures>`*`,`*`<numeric>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`infIsNA(`*`<QFeatures>`*`,`*`<character>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`nNA(`*`<SummarizedExperiment>`*`,`*`<missing>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`nNA(`*`<QFeatures>`*`,`*`<integer>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`nNA(`*`<QFeatures>`*`,`*`<numeric>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`nNA(`*`<QFeatures>`*`,`*`<character>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`filterNA(`*`<SummarizedExperiment>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  [`filterNA(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-missing-data.md)
  : Managing missing data
- [`impute`](https://rformassspectrometry.github.io/QFeatures/reference/impute.md)
  [`impute(`*`<SummarizedExperiment>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/impute.md)
  [`impute(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/impute.md)
  : Quantitative proteomics data imputation

## QFeatures processing

Quantitative data processing

- [`logTransform(`*`<SummarizedExperiment>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md)
  [`logTransform(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md)
  [`scaleTransform(`*`<SummarizedExperiment>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md)
  [`scaleTransform(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md)
  [`normalize(`*`<SummarizedExperiment>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md)
  [`normalize(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md)
  [`sweep(`*`<SummarizedExperiment>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md)
  [`sweep(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-processing.md)
  : QFeatures processing

## Feature aggregations

Feature aggregations

- [`aggregateFeatures(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md)
  [`aggregateFeatures(`*`<SummarizedExperiment>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md)
  [`adjacencyMatrix(`*`<QFeatures>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md)
  [`` `adjacencyMatrix<-`() ``](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md)
  [`aggcounts(`*`<SummarizedExperiment>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-aggregate.md)
  : Aggregate assays' quantitative features

## Filtering

Feature filtering

- [`VariableFilter()`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-filtering.md)
  [`filterFeatures(`*`<QFeatures>`*`,`*`<AnnotationFilter>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-filtering.md)
  [`filterFeatures(`*`<QFeatures>`*`,`*`<formula>`*`)`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-filtering.md)
  : Filter features based on their rowData
- [`subsetByFeature`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-subsetBy.md)
  [`subsetByFeature,QFeatures,character-method`](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-subsetBy.md)
  : Subset by feature name

## Example data

Example data

- [`feat1`](https://rformassspectrometry.github.io/QFeatures/reference/feat1.md)
  : Feature example data

- [`feat3`](https://rformassspectrometry.github.io/QFeatures/reference/feat3.md)
  :

  Example `QFeatures` object after processing

- [`hlpsms`](https://rformassspectrometry.github.io/QFeatures/reference/hlpsms.md)
  : hyperLOPIT PSM-level expression data

## DataFrame

DataFrame functions

- [`reduceDataFrame()`](https://rformassspectrometry.github.io/QFeatures/reference/reduceDataFrame.md)
  [`expandDataFrame()`](https://rformassspectrometry.github.io/QFeatures/reference/reduceDataFrame.md)
  :

  Reduces and expands a `DataFrame`

- [`unfoldDataFrame()`](https://rformassspectrometry.github.io/QFeatures/reference/unfoldDataFrame.md)
  : Unfold a data frame

## Generic functions

Package generics

- [`AllGenerics`](https://rformassspectrometry.github.io/QFeatures/reference/AllGenerics.md)
  : Placeholder for generics functions documentation
