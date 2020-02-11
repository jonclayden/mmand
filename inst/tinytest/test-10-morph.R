source("expectation.R")

library(loder)
fan <- readPng(system.file("images", "fan-small.png", package="mmand"))

# Binary mathematical morphology
data <- c(0,0,1,0,0,0,1,1,1,0,0)
kernel <- c(1,1,1)
expect_true(binary(data))
expect_equal(erode(data,kernel), c(0,0,0,0,0,0,0,1,0,0,0))
expect_equal(dilate(data,kernel), c(0,1,1,1,0,1,1,1,1,1,0))
expect_equal(skeletonise(data,kernel), c(0,0,1,0,0,0,0,1,0,0,0))

# Odd kernels: asymmetric and zero-origin
expect_equal(erode(data,c(0,1,1)), c(0,0,0,0,0,0,1,1,0,0,0))
expect_equal(dilate(data,c(0,1,1)), c(0,0,1,1,0,0,1,1,1,1,0))
expect_equal(erode(data,c(1,0,1)), c(0,0,0,0,0,0,0,1,0,0,0))
expect_equal(dilate(data,c(1,0,1)), c(0,1,0,1,0,1,1,1,1,1,0))

data <- matrix(0, nrow=3, ncol=3)
data[2,2] <- 1
kernel <- shapeKernel(c(3,3), type="diamond")
expect_equal(neighbourhood(data,3)$offsets, -4:4)
expect_equal_to_reference(dilate(data,kernel), "dilate_2d_bin.rds")


# Greyscale mathematical morphology
data <- c(0,0,0.5,0,0,0,0.2,0.5,0.3,0,0)
kernel <- c(1,1,1)
expect_false(binary(data))
expect_equal(erode(data,kernel), c(0,0,0,0,0,0,0,0.2,0,0,0))
expect_equal(dilate(data,kernel), c(0,0.5,0.5,0.5,0,0.2,0.5,0.5,0.5,0.3,0))

kernel <- c(0.5,1,0.5)
expect_equal(erode(data,kernel), c(-1,-1,-0.5,-1,-1,-1,-0.8,-0.5,-0.7,-1,-1))
expect_equal(dilate(data,kernel), c(1,1,1.5,1,1,1,1.2,1.5,1.3,1,1))

kernel <- shapeKernel(c(3,3), type="diamond")
expect_equal_to_reference(erode(fan,kernel), "fan_eroded.rds")
expect_equal_to_reference(dilate(fan,kernel), "fan_dilated.rds")
expect_equal_to_reference(closing(fan,kernel), "fan_opened.rds")
expect_equal_to_reference(opening(fan,kernel), "fan_closed.rds")


# Smoothing and filtering
data <- matrix(0, nrow=3, ncol=3)
data[2,2] <- 1
expect_equal_to_reference(gaussianSmooth(data,c(1,1)), "2d_smooth_small.rds")

data <- matrix(0, nrow=7, ncol=7)
data[4,4] <- 1
expect_equal_to_reference(gaussianSmooth(data,c(1,1)), "2d_smooth_large.rds")

kernel <- shapeKernel(c(3,3), type="diamond")
expect_equal_to_reference(meanFilter(fan,kernel), "fan_mean_filtered.rds")
expect_equal_to_reference(medianFilter(fan,kernel), "fan_median_filtered.rds")

expect_equal_to_reference(sobelFilter(fan), "fan_sobel_filtered.rds")


# Binarising and thresholding
data <- c(0.1, 0.05, 0.95, 0.85, 0.15, 0.9)
expect_equal(binarise(data)[3], 1)
expect_equal(threshold(data,0.5)[3], 1)
expect_equal(threshold(data,0.5,binarise=FALSE)[3], 0.95)
expect_equal(threshold(data,method="kmeans")[3], 1)
