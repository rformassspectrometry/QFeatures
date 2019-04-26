## Taken from SummarizedExperiment to reproduce the
## SummarizedExperiment,show method - need to request to put it out of
## the method's body for re-use.
scat <- function(fmt, vals = character(), exdent = 2, ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse = " ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent=exdent, ...), sep = "\n")
}
