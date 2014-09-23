library(png)
lena <- readPNG(system.file("images", "lena-small.png", package="mmand"))

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
    expect_that(erode(lena,kernel), equals_reference("lena_eroded.rds"))
    expect_that(dilate(lena,kernel), equals_reference("lena_dilated.rds"))
    expect_that(closing(lena,kernel), equals_reference("lena_opened.rds"))
    expect_that(opening(lena,kernel), equals_reference("lena_closed.rds"))
})

test_that("smoothing and filtering operations work", {
    data <- matrix(0, nrow=3, ncol=3)
    data[2,2] <- 1
    expect_that(gaussianSmooth(data,c(1,1)), equals_reference("2d_smooth_small.rds"))
    
    data <- matrix(0, nrow=7, ncol=7)
    data[4,4] <- 1
    expect_that(gaussianSmooth(data,c(1,1)), equals_reference("2d_smooth_large.rds"))
    
    kernel <- shapeKernel(c(3,3), type="diamond")
    expect_that(meanFilter(lena,kernel), equals_reference("lena_mean_filtered.rds"))
    expect_that(medianFilter(lena,kernel), equals_reference("lena_median_filtered.rds"))
})
