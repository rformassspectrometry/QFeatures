# One or two levels in Features

One or two assay levels could be considered in Features:

- one level: each SE contains only a single assay, and when an SE is
  processed (log-transformed, normalised, ...) in a way that its
  dimensions stay the same, a new SE is created and added to the
  Features object.
  
- two level: SEs can contain multiple assays, and when an SE is
  processed (log-transformed, normalised, ...) in a way that its
  dimensions stay the same, a new assay is added to that SE.
  
This
[question](https://stat.ethz.ch/pipermail/bioc-devel/2020-January/016096.html)
on the bioc-devel list ask for advice on SE processing, and whether a
new SE or new assay in the original SE should be preferred. While the
letter is arguably more elegant, and is also used in
SingleAssayExperiment pipelines, it doesn't seem to be the case when
using SummarizedExperiments.

As for features (or [MultiAssayExperiments in
general](https://github.com/waldronlab/MultiAssayExperiment/issues/266)),
the two-level approach isn't readily available out-of-the-box, and
would require additional developments:

- Every functoin that operates on an SE of a Features object would
  need to allow the user to specify which assay to use (and/or by
  default use the latest one).
  
- The `show,Features` method would need to display the number/names of
  the assays in each SE to make these two levels explicit.

Despite the elegant of the two-level option, it seems that the
additional development isn't warranted at this time.


The [`updateAssay`
function](https://github.com/rformassspectrometry/Features/issues/37)
was originally intended for the two-level approach, i.e. to add an
assay to an SE. This is not considered anymore (for now, at least).
