library("Features")

psms <- matrix(1:20, ncol = 2)
colnames(psms) <- paste0("S", 1:2)
rowdata <- DataFrame(Sequence = c("SYGFNAAR", "SYGFNAAR", "SYGFNAAR", "ELGNDAYK",
                                  "ELGNDAYK", "ELGNDAYK", "IAEESNFPFIK",
                                  "IAEESNFPFIK", "IAEESNFPFIK", "IAEESNFPFIK"),
                     Protein = c("ProtA", "ProtA", "ProtA", "ProtA", "ProtA",
                                 "ProtA", "ProtB", "ProtB", "ProtB", "ProtB"),
                     Var = 1:10)
rownames(rowdata) <- rownames(psms) <- paste0("PSM", 1:10)
coldata <- DataFrame(Group = 1:2)
rownames(coldata) <- colnames(psms)


psms <- SummarizedExperiment(psms, rowData = rowdata)

feat1 <- Features(list(psmsm = psms), colData = coldata)

save(feat1, file = "../../data/feat1.rda", compress = "xz", compression_level = 9)
