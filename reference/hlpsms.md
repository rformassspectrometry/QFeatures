# hyperLOPIT PSM-level expression data

A `data.frame` with PSM-level quantitation data by Christoforou *et al.*
(2016). This is the first replicate of a spatial proteomics dataset from
a hyperLOPIT experimental design on Mouse E14TG2a embryonic stem cells.
Normalised intensities for proteins for TMT 10-plex labelled fractions
are available for 3 replicates acquired in MS3 mode using an Orbitrap
Fusion mass-spectrometer.

The variable names are

- X126, X127C, X127N, X128C, X128N, X129C, X129N, X130C, X130N and X131:
  the 10 TMT tags used to quantify the peptides along the density
  gradient.

- Sequence: the peptide sequence.

- ProteinDescriptions: the description of the protein this peptide was
  associated to.

- NbProteins: the number of proteins in the protein group.

- ProteinGroupAccessions: the main protein accession number in the
  protein group.

- Modifications: post-translational modifications identified in the
  peptide.

- qValue: the PSM identification q-value.

- PEP: the PSM posterior error probability.

- IonScore: the Mascot ion identification score.

- NbMissedCleavages: the number of missed cleavages in the peptide.

- IsolationInterference: the calculated precursor ion isolation
  interference.

- IonInjectTimems: the ions injection time in milli-seconds.

- Intensity: the precursor ion intensity.

- Charge: the peptide charge.

- mzDa: the peptide mass to charge ratio, in Daltons.

- MHDa: the peptide mass, in Daltons.

- DeltaMassPPM: the difference in measure and calculated mass, in parts
  per millions.

- RTmin: the peptide retention time, in minutes.

- markers: localisation for well known sub-cellular markers. QFeatures
  of unknown location are encode as `"unknown"`.

For further details, install the `pRolocdata` package and see
`?hyperLOPIT2015`.

## Usage

``` r
hlpsms
```

## Format

An object of class `data.frame` with 3010 rows and 28 columns.

## Source

The `pRolocdata` package: <http://bioconductor.org/packages/pRolocdata/>

## References

*A draft map of the mouse pluripotent stem cell spatial proteome*
Christoforou A, Mulvey CM, Breckels LM, Geladaki A, Hurrell T, Hayward
PC, Naake T, Gatto L, Viner R, Martinez Arias A, Lilley KS. Nat Commun.
2016 Jan 12;7:8992. doi: 10.1038/ncomms9992. PubMed PMID: 26754106;
PubMed Central PMCID: PMC4729960.

## See also

See
[QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.md)
to import this data using the
[`readQFeatures()`](https://rformassspectrometry.github.io/QFeatures/reference/readQFeatures.md)
function.
