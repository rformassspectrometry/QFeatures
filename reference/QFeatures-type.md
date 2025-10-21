# Set and Get QFeatures Type

Developer-level functions to set and retrieve the type of a `QFeatures`
object. This type can help internal methods adapt their behaviour to the
structure of the data.

## Usage

``` r
setQFeaturesType(object, type)

getQFeaturesType(object)

validQFeaturesTypes()
```

## Arguments

- object:

  An instance of class
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md).

- type:

  `character(1)` defining the type of the QFeatures. Must be one of the
  values returned by `validQFeaturesTypes()`.

## Value

- `setQFeaturesType()`: returns the updated `QFeatures` object with the
  type stored in its metadata.

- `getQFeaturesType()`: returns a character string indicating the type
  of the `QFeatures` object. If no type is explicitly set, it is
  inferred from the class of the experiments. If the QFeatures contains
  any `SingleCellExperiment` objects, the type is set to "scp".
  Otherwise, it is set to "bulk".

- `validQFeaturesTypes()`: character vector of valid QFeatures types.

## Details

These functions control an internal metadata slot (`._type`) used to
distinguish between different structural uses of `QFeatures` objects.
This slot is directly accessible with `metadata(object)[["._type"]]`.

## Note

The `QFeatures` type slot was introduced because, in the context of the
`scp` package, we found that `SingleCellExperiment` objects were slower
than `SummarizedExperiment` objects ([GH issue:
scp#83](https://github.com/UCLouvain-CBIO/scp/issues/83)). As a result,
we started using `SummarizedExperiment` objects within `scp`. However,
to retain information about the type of data being handled, we
introduced the `QFeatures` type slot. This slot is, for example, used in
the `show` method of `QFeatures`.

## Warning

These functions are intended for package developers and internal use.
End users should typically not call them directly.
