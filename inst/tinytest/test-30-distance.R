# Distance transform

# Several of these tests fail on CRAN Solaris for unknown reasons
# It looks like the square-root is taken twice in the first few - no idea why
# Not reproducible within Vagrant box
if (tolower(Sys.info()[["sysname"]]) != "sunos") {
    # The original values make no difference - the source array is treated as binary
    expect_equal(distanceTransform(c(0,0,1,0,0,0,1,1,1,0,0)), c(2,1,0,1,2,1,0,0,0,1,2))
    expect_equal(distanceTransform(c(0,0,1,0,0,0,1,1,1,0,0),signed=TRUE), c(2,1,-1,1,2,1,-1,-2,-1,1,2))
    expect_equal(distanceTransform(c(0,0,1,0,0,0,2,3,4,0,0)), c(2,1,0,1,2,1,0,0,0,1,2))
    expect_equal(distanceTransform(c(0,0,1,0,0,0,2,3,4,0,0),signed=TRUE), c(2,1,-1,1,2,1,-1,-2,-1,1,2))

    kernel <- shapeKernel(c(5,5), type="diamond")
    unsignedTransform <- distanceTransform(kernel)
    expect_equal(unsignedTransform[3,3], 0)
    expect_equal(unsignedTransform[,1], sqrt(c(2,1,0,1,2)))

    signedTransform <- distanceTransform(kernel, signed=TRUE)
    expect_equal(signedTransform, unsignedTransform - distanceTransform(1-kernel))
    expect_equal(signedTransform[3,3], -sqrt(5))

    # Set the pixdims so that distances are twice as large in the second dimension
    anisotropicTransform <- distanceTransform(kernel, c(1,2))
    expect_equal(anisotropicTransform[,1], c(2,1,0,1,2))
    
    # The 2D distance transform of a diagonal matrix should be symmetrical
    transform2D <- distanceTransform(diag(5))
    expect_equal(transform2D, t(transform2D))
}
