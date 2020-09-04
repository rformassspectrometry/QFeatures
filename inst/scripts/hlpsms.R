library("pRolocdata")
library("MSnbase")

data("hyperLOPIT2015ms3r1psm")

vars <- c("Sequence", "Protein.Descriptions", "X..Proteins",
          "Protein.Group.Accessions", "Modifications", "q.Value", "PEP",
          "IonScore", "X..Missed.Cleavages", "Isolation.Interference....",
          "Ion.Inject.Time..ms.", "Intensity", "Charge", "m.z..Da.", "MH...Da.",
          "Delta.Mass..PPM.", "RT..min.", "markers")

hlpsms  <- ms2df(selectFeatureData(hyperLOPIT2015ms3r1psm, fcol = vars))

names(hlpsms) <- gsub("\\.", "", names(hlpsms))
names(hlpsms) <- gsub("XProteins", "NbProteins", names(hlpsms))
names(hlpsms) <- gsub("XMissedCleavages", "NbMissedCleavages", names(hlpsms))

for (i in seq_along(hlpsms))
    if (is(hlpsms[[i]], "factor"))
        hlpsms[[i]] <- as.character(hlpsms[[i]])

## 2020-08-12: subset the data.frame to save space, but keep the STAT1 and STAT3
## features (used in the vignette)
stats <- grep("STAT", hlpsms$ProteinDescriptions)
set.seed(123)
k <- sample(nrow(hlpsms), 3000)
hlpsms <- hlpsms[sort(union(k, stats)), ]

save(hlpsms, file = "../../data/hlpsms.rda", compress = "xz", compression_level = 9)
