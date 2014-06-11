source("matches_file.R")

context("Generating and sampling kernels")

test_that("standard kernel arrays can be created", {
    expect_that(shapeKernel(3), matches_file("line_kernel.rds"))
    expect_that(shapeKernel(c(5,5),type="box"), matches_file("box_kernel.rds"))
    expect_that(shapeKernel(c(5,5),type="disc"), matches_file("disc_kernel.rds"))
    expect_that(shapeKernel(c(5,5),type="diamond"), matches_file("diamond_kernel.rds"))
    
    expect_that(gaussianKernel(0.5), matches_file("gaussian_kernel_1d.rds"))
    expect_that(gaussianKernel(c(0.5,0.5)), matches_file("gaussian_kernel_2d.rds"))
    expect_that(gaussianKernel(c(0.5,0.3)), matches_file("gaussian_kernel_2d_anis.rds"))
})

test_that("standard kernel functions can be created", {
    expect_that(boxKernel(), matches_file("box_function.rds"))
    expect_that(triangleKernel(), matches_file("triangle_function.rds"))
    expect_that(mitchellNetravaliKernel(), matches_file("mn_function.rds"))
})

test_that("we can sample from kernel functions", {
    expect_that(sampleKernelFunction(boxKernel(),seq(-1,1,0.5)), equals(c(0,1,1,1,0)))
    expect_that(sampleKernelFunction(triangleKernel(),seq(-1,1,0.5)), equals(c(0,0.5,1,0.5,0)))
    expect_that(sampleKernelFunction(mitchellNetravaliKernel(0,1),seq(-1,1,0.5)), equals(c(0,0.625,1,0.625,0)))
})
