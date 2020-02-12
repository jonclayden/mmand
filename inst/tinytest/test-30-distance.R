# Distance transform
expect_equal(distanceTransform(c(0,0,1,0,0,0,1,1,1,0,0)), c(2,1,0,1,2,1,0,0,0,1,2))

kernel <- shapeKernel(c(5,5), type="diamond")
unsignedTransform <- distanceTransform(kernel)
expect_equal(unsignedTransform[3,3], 0)
expect_equal(unsignedTransform[,1], sqrt(c(2,1,0,1,2)))

signedTransform <- unsignedTransform - distanceTransform(1-kernel)
expect_equal(signedTransform[3,3], -sqrt(5))
