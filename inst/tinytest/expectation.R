expect_equal_to_reference <- function (expr, file, ...)
{
    if (file.exists(file))
        tinytest::expect_equivalent(expr, readRDS(file), ...)
    else
        tinytest::expect_null(saveRDS(expr, file))
}
