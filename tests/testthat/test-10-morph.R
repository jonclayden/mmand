library(png)
fan <- readPNG(system.file("images", "fan-small.png", package="mmand"))

context("Mathematical morphology and filtering")

test_that("binary mathematical morphology operations work", {
    data <- c(0,0,1,0,0,0,1,1,1,0,0)
    kernel <- c(1,1,1)
    expect_that(binary(data), is_true())
    expect_that(erode(data,kernel), equals(c(0,0,0,0,0,0,0,1,0,0,0)))
    expect_that(dilate(data,kernel), equals(c(0,1,1,1,0,1,1,1,1,1,0)))
    
    data <- matrix(0, nrow=3, ncol=3)
    data[2,2] <- 1
    kernel <- shapeKernel(c(3,3), type="diamond")
    expect_that(dilate(data,kernel), equals_reference("dilate_2d_bin.rds"))
})

test_that("greyscale mathematical morphology operations work", {
    data <- c(0,0,0.5,0,0,0,0.2,0.5,0.3,0,0)
    kernel <- c(1,1,1)
    expect_that(binary(data), is_false())
    expect_that(erode(data,kernel), equals(c(0,0,0,0,0,0,0,0.2,0,0,0)))
    expect_that(dilate(data,kernel), equals(c(0,0.5,0.5,0.5,0,0.2,0.5,0.5,0.5,0.3,0)))
    
    kernel <- c(0.5,1,0.5)
    expect_that(erode(data,kernel), equals(c(-1,-1,-0.5,-1,-1,-1,-0.8,-0.5,-0.7,-1,-1)))
    expect_that(dilate(data,kernel), equals(c(1,1,1.5,1,1,1,1.2,1.5,1.3,1,1)))
    
    kernel <- shapeKernel(c(3,3), type="diamond")
    expect_that(erode(fan,kernel), equals_reference("fan_eroded.rds"))
    expect_that(dilate(fan,kernel), equals_reference("fan_dilated.rds"))
    expect_that(closing(fan,kernel), equals_reference("fan_opened.rds"))
    expect_that(opening(fan,kernel), equals_reference("fan_closed.rds"))
})

test_that("smoothing and filtering operations work", {
    data <- matrix(0, nrow=3, ncol=3)
    data[2,2] <- 1
    expect_that(gaussianSmooth(data,c(1,1)), equals_reference("2d_smooth_small.rds"))
    
    data <- matrix(0, nrow=7, ncol=7)
    data[4,4] <- 1
    expect_that(gaussianSmooth(data,c(1,1)), equals_reference("2d_smooth_large.rds"))
    
    kernel <- shapeKernel(c(3,3), type="diamond")
    expect_that(meanFilter(fan,kernel), equals_reference("fan_mean_filtered.rds"))
    expect_that(medianFilter(fan,kernel), equals_reference("fan_median_filtered.rds"))
})
