# Links between Assays

Links between assays within a
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
object are handled by an `AssayLinks` object. It is composed by a list
of `AssayLink` instances.

## Usage

``` r
# S4 method for class 'AssayLink'
show(object)

# S4 method for class 'AssayLinks'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'AssayLink'
updateObject(object, ..., verbose = FALSE)

AssayLink(name, from = NA_character_, fcol = NA_character_, hits = Hits())

AssayLinks(..., names = NULL)

assayLink(x, i)

assayLinks(x, i)

# S4 method for class 'AssayLink,character,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'AssayLinks,list,ANY,ANY'
x[i, j, ..., drop = TRUE]

addAssayLink(object, from, to, varFrom, varTo)

addAssayLinkOneToOne(object, from, to)
```

## Arguments

- object:

  An `AssayLink` object to show.

- ...:

  A set of `AssayLink` objects or a list thereof.

- verbose:

  logical (default FALSE) whether to print extra messages

- name:

  A mandatory name of the assay(s).

- from:

  A [`character()`](https://rdrr.io/r/base/character.html) or
  [`integer()`](https://rdrr.io/r/base/integer.html) indicating which
  assay(s) to link from in `object`

- fcol:

  The feature variable of the parent assay used to generate the current
  assay (used in `aggregateFeatures`). `NA_character_`, if not
  applicable.

- hits:

  An object of class
  [S4Vectors::Hits](https://rdrr.io/pkg/S4Vectors/man/Hits-class.html)
  matching the features of two assays.

- names:

  A [`character()`](https://rdrr.io/r/base/character.html) of
  `AssayLink` names. If provided, `...` are ignored, and `names` is used
  to create an `AssayLinks` object with `AssayLink` instances with names
  `names`.

- x:

  An instance of class
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md).

- i:

  The index or name of the assay whose `AssayLink` and parents
  `AssayLink` instances are to be returned. For `[`, the feature names
  to filter on.

- j:

  ignored.

- drop:

  ignored.

- to:

  A `character(1)` or `integer(1)` indicating which assay to link to in
  `object`

- varFrom:

  A [`character()`](https://rdrr.io/r/base/character.html) indicating
  the feature variable(s) to use to match the `from` assay(s) to the
  `to` assay. `varFrom` must have the same length as `from` and is
  assumed to be ordered as `from`.

- varTo:

  A `character(1)` indicating the feature variable to use to match the
  `to` assay to the `from` assay(s).

## Value

`assayLink` returns an instance of class `AssayLink`.

`assayLinks` returns an instance of class `AssayLinks`.

## Constructors

Object can be created with the `AssayLink()` and `AssayLinks()`
constructors.

## Methods and functions

- `assayLink(x, i)` accesses the AssayLink at position `i` or with name
  `i` in the
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  object `x`.

- `parentAssayLinks(x, i, recursive = FALSE)` accesses the parent(s)
  `AssayLinks` or assay with index or name `i`.

## Creating links between assays

- `addAssayLink` takes a parent assay and a child assay contained in the
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  object and creates a link given a matching feature variable in each
  assay's `rowData`. `addAssayLink` also allows to link an assay from
  multiple parent assays (see examples below).

- `addAssayLinkOneToOne` links two assays contained in the
  [QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
  object. The parent assay and the child assay must have the same size
  and contain the same rownames (a different ordering is allowed). The
  matching is performed based on the row names of the assays, instead of
  a supplied variable name in `rowData`. Providing multiple parents is
  not supported.

## Examples

``` r

##-----------------------------
## Creating an AssayLink object
##-----------------------------

al1 <- AssayLink(name = "assay1")
al1
#> AssayLink for assay <assay1>
#> [from:NA|fcol:NA|hits:0]

##------------------------------
## Creating an AssayLinks object
##------------------------------

AssayLinks(al1)
#> AssayLinks of length 1
#> names(1): assay1

al2 <- AssayLinks(names = c("Assay1", "Assay2"))
al2
#> AssayLinks of length 2
#> names(2): Assay1 Assay2

##---------------------------------------
## Adding an AssayLink between two assays
##---------------------------------------

## create a QFeatures object with 2 (identical) assays
## see also '?QFeatures'
se <- SummarizedExperiment(matrix(runif(20), ncol = 2,
                                  dimnames = list(LETTERS[1:10],
                                                  letters[1:2])),
                           rowData = DataFrame(ID = 1:10))
ft <- QFeatures(list(assay1 = se, assay2 = se))

## assay1 and assay2 are not linked
assayLink(ft, "assay2") ## 'from' is NA
#> AssayLink for assay <assay2>
#> [from:NA|fcol:NA|hits:0]
assayLink(ft, "assay1") ## 'from' is NA
#> AssayLink for assay <assay1>
#> [from:NA|fcol:NA|hits:0]

## Suppose assay2 was generated from assay1 and the feature variable
## 'ID' keeps track of the relationship between the two assays
ftLinked <- addAssayLink(ft, from = "assay1", to = "assay2",
                         varFrom = "ID", varTo = "ID")
assayLink(ftLinked, "assay2")
#> AssayLink for assay <assay2>
#> [from:assay1|fcol:ID|hits:10]

## For one-to-one relationships, you can also use
ftLinked <- addAssayLinkOneToOne(ft, from = "assay1", to = "assay2")
assayLink(ftLinked, "assay2")
#> AssayLink for assay <assay2>
#> [from:assay1|fcol:._oneToOne|hits:10]

##----------------------------------------
## Adding an AssayLink between more assays
##----------------------------------------

## An assay can also be linked to multiple parent assays
## Create a QFeatures object with 2 parent assays and 1 child assay
ft <- QFeatures(list(parent1 = se[1:6, ], parent2 = se[4:10, ], child = se))
ft <- addAssayLink(ft, from = c("parent1", "parent2"), to = "child",
                   varFrom = c("ID", "ID"), varTo = "ID")
assayLink(ft, "child")
#> AssayLink for assay <child>
#> [from:parent1,parent2|fcol:ID,ID|hits:6,7]
```
