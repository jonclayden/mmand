matches_file <- function (file, compare = equals, label = NULL, ...)
{
    if (file.exists(file)) {
        reference <- readRDS(file)
        compare <- match.fun(compare)
        if (is.null(label))
            label <- paste("reference from", file)
        compare(reference, label=label, ...)
    } else {
        return (function(actual) {
            saveRDS(actual, file)
            expectation(TRUE, "should never fail", "saved to file")
        })
    }
}
