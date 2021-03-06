fan <- loder::readPng(system.file("images", "fan-small.png", package="mmand"))

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
expect_equal(dilate(data,kernel), readRDS("dilate_2d_bin.rds"))

# Different skeletonisation methods (in 2D)
data <- shapeKernel(c(5,5), type="diamond")
kernel <- shapeKernel(c(3,3))
expect_equal(skeletonise(data,kernel,method="lantuejoul")[,3], c(1,0,1,0,1))
expect_equal(skeletonise(data,kernel,method="beucher")[,3], c(1,1,1,1,1))
expect_equal(skeletonise(data,kernel,method="hitormiss")[,3], c(0,0,1,1,1))

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
expect_equal(erode(fan,kernel), readRDS("fan_eroded.rds"))
expect_equal(dilate(fan,kernel), readRDS("fan_dilated.rds"))
expect_equal(closing(fan,kernel), readRDS("fan_opened.rds"))
expect_equal(opening(fan,kernel), readRDS("fan_closed.rds"))


# Smoothing and filtering
data <- matrix(0, nrow=3, ncol=3)
data[2,2] <- 1
expect_equal(gaussianSmooth(data,c(1,1)), readRDS("2d_smooth_small.rds"))

data <- matrix(0, nrow=7, ncol=7)
data[4,4] <- 1
expect_equal(gaussianSmooth(data,c(1,1)), readRDS("2d_smooth_large.rds"))

kernel <- shapeKernel(c(3,3), type="diamond")
expect_equal(meanFilter(fan,kernel), readRDS("fan_mean_filtered.rds"))
expect_equal(medianFilter(fan,kernel), readRDS("fan_median_filtered.rds"))

expect_equal(sobelFilter(fan), readRDS("fan_sobel_filtered.rds"))


# Binarising and thresholding
data <- c(0.1, 0.05, 0.95, 0.85, 0.15, 0.9)
expect_equal(binarise(data)[3], 1)
expect_equal(threshold(data,0.5)[3], 1)
expect_equal(threshold(data,0.5,binarise=FALSE)[3], 0.95)
expect_equal(threshold(data,method="kmeans")[3], 1)
