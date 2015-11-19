context("Kernel generation and sampling")

test_that("standard kernel arrays can be created", {
    expect_that(shapeKernel(3), equals_reference("line_kernel.rds"))
    expect_that(shapeKernel(c(5,5),type="box"), equals_reference("box_kernel.rds"))
    expect_that(shapeKernel(c(5,5),type="disc"), equals_reference("disc_kernel.rds"))
    expect_that(shapeKernel(c(5,5),type="diamond"), equals_reference("diamond_kernel.rds"))
    
    expect_that(gaussianKernel(0.5), equals_reference("gaussian_kernel_1d.rds"))
    expect_that(gaussianKernel(c(0.5,0.5)), equals_reference("gaussian_kernel_2d.rds"))
    expect_that(gaussianKernel(c(0.5,0.5),normalised=FALSE), equals_reference("gaussian_kernel_2d_unnorm.rds"))
    expect_that(gaussianKernel(c(0.5,0.3)), equals_reference("gaussian_kernel_2d_anis.rds"))
    
    expect_that(sobelKernel(1), equals(kernelArray(c(1,0,-1))))
    expect_that(sobelKernel(1,0), equals(kernelArray(c(1,2,1)/4)))
    expect_that(sobelKernel(2), equals_reference("sobel_kernel_2d.rds"))
})

test_that("standard kernel functions can be created", {
    expect_that(boxKernel(), equals_reference("box_function.rds"))
    expect_that(triangleKernel(), equals_reference("triangle_function.rds"))
    expect_that(mitchellNetravaliKernel(), equals_reference("mn_function.rds"))
})

test_that("we can sample from kernel functions", {
    expect_that(sampleKernelFunction(boxKernel(),seq(-1,1,0.5)), equals(c(0,1,1,1,0)))
    expect_that(sampleKernelFunction(triangleKernel(),seq(-1,1,0.5)), equals(c(0,0.5,1,0.5,0)))
    expect_that(sampleKernelFunction(mitchellNetravaliKernel(0,1),seq(-1,1,0.5)), equals(c(0,0.625,1,0.625,0)))
})
