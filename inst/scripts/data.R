library("pRolocdata")

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

save(hlpsms, file = "../../data/hlpsms.rda", compress = "xz", compression_level = 9)