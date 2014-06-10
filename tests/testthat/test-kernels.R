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

test_that("standard array kernels can be created", {
    expect_that(shapeKernel(3), matches_file("line_kernel.rds"))
    expect_that(shapeKernel(c(5,5),type="box"), matches_file("box_kernel.rds"))
    expect_that(shapeKernel(c(5,5),type="diamond"), matches_file("disc_kernel.rds"))
    expect_that(shapeKernel(c(5,5),type="diamond"), matches_file("diamond_kernel.rds"))
})
