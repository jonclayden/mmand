# Standard kernel arrays
expect_equal_to_reference(shapeKernel(3), "line_kernel.rds")
expect_equal_to_reference(shapeKernel(c(5,5),type="box"), "box_kernel.rds")
expect_equal_to_reference(shapeKernel(c(5,5),type="disc"), "disc_kernel.rds")
expect_equal_to_reference(shapeKernel(c(5,5),type="diamond"), "diamond_kernel.rds")

expect_equal_to_reference(gaussianKernel(0.5), "gaussian_kernel_1d.rds")
expect_equal_to_reference(gaussianKernel(c(0.5,0.5)), "gaussian_kernel_2d.rds")
expect_equal_to_reference(gaussianKernel(c(0.5,0.5),normalised=FALSE), "gaussian_kernel_2d_unnorm.rds")
expect_equal_to_reference(gaussianKernel(c(0.5,0.3)), "gaussian_kernel_2d_anis.rds")

expect_equal(sobelKernel(1), kernelArray(c(1,0,-1)))
expect_equal(sobelKernel(1,0), kernelArray(c(1,2,1)/4))
expect_equal_to_reference(sobelKernel(2), "sobel_kernel_2d.rds")


# Standard kernel functions
expect_equal_to_reference(boxKernel(), "box_function.rds")
expect_equal_to_reference(triangleKernel(), "triangle_function.rds")
expect_equal_to_reference(mitchellNetravaliKernel(), "mn_function.rds")


# Sampling from kernel functions
expect_equal(sampleKernelFunction(boxKernel(),seq(-1,1,0.5)), c(0,1,1,1,0))
expect_equal(sampleKernelFunction(triangleKernel(),seq(-1,1,0.5)), c(0,0.5,1,0.5,0))
expect_equal(sampleKernelFunction(mitchellNetravaliKernel(0,1),seq(-1,1,0.5)), c(0,0.625,1,0.625,0))


# Type testing
expect_true(isKernel(boxKernel()))
expect_true(isKernelFunction(boxKernel()))
expect_false(isKernelArray(boxKernel()))
expect_true(isKernel(shapeKernel(3)))
expect_false(isKernelFunction(shapeKernel(3)))
expect_true(isKernelArray(shapeKernel(3)))


# Plotting
plot(boxKernel())
plot(shapeKernel(3))
