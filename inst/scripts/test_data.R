library("Features")

set.seet(1)
## Creating psms as <double> rather than <integer>
psms <- matrix(seq(1, 20, by = 1.0), ncol = 2)
colnames(psms) <- paste0("S", 1:2)
rowdata <- DataFrame(Sequence = c("SYGFNAAR", "SYGFNAAR", "SYGFNAAR", "ELGNDAYK",
                                  "ELGNDAYK", "ELGNDAYK", "IAEESNFPFIK",
                                  "IAEESNFPFIK", "IAEESNFPFIK", "IAEESNFPFIK"),
                     Protein = c("ProtA", "ProtA", "ProtA", "ProtA", "ProtA",
                                 "ProtA", "ProtB", "ProtB", "ProtB", "ProtB"),
                     Var = 1:10,
                     location = c(rep("Mitochondrion", 6), rep("unknown", 4)),
                     pval = round(runif(10, max = 0.1, min = 0), 3))
rownames(rowdata) <- rownames(psms) <- paste0("PSM", 1:10)
coldata <- DataFrame(Group = 1:2)
rownames(coldata) <- colnames(psms)


psms <- SummarizedExperiment(psms, rowData = rowdata)

feat1 <- Features(list(psms = psms), colData = coldata)

save(feat1, file = "../../data/feat1.rda", compress = "xz", compression_level = 9)
