B <- loder::readPng(system.file("images", "B.png", package="mmand"))

# We need to avoid creating plot files within installed package directories
setwd(tempdir())
expect_null(display(B))

output <- capture.output(sketch(B))
expect_true(any(grepl("@@@@@", output, fixed=TRUE)))
expect_false(all(grepl("@@@@@", output, fixed=TRUE)))
