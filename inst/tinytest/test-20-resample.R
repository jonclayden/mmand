# Image resampling and rescaling
expect_equal(resample(c(0,0,1,0,0),seq(0.75,5.25,0.5),boxKernel()), c(0,0,0,0,1,1,0,0,0,0))
expect_equal(resample(c(0,0,1,0,0),seq(0.75,5.25,0.5),triangleKernel()), c(0,0,0,0.25,0.75,0.75,0.25,0,0,0))

expect_equal(resample(1:5,1:5,mitchellNetravaliKernel()), as.numeric(1:5))
expect_equal(resample(1:5,3.5,mitchellNetravaliKernel()), 3.5)

data <- matrix(1:9, nrow=3, ncol=3)
point <- matrix(c(1.5,2.5), nrow=1)
expect_equal(resample(data,point,boxKernel()), 8)
expect_equal(resample(data,point,triangleKernel()), 6)
expect_equal(resample(data,point,mitchellNetravaliKernel()), 6)
expect_equal(resample(data,point,lanczosKernel()), 5.527249,tol=0.001)

points <- point %x% matrix(1,4,1)
expect_equal(resample(data,points,mitchellNetravaliKernel()), c(6,6,6,6))

grid <- list(c(1.5,2.5),c(1.5,2.5))
expect_equal(resample(data,grid,boxKernel()), matrix(c(5,6,8,9),2,2))
expect_equal(resample(data,grid,triangleKernel()), matrix(c(3,4,6,7),2,2))
expect_equal(resample(data,grid,mitchellNetravaliKernel()), matrix(c(3,4,6,7),2,2))

expect_equal(rescale(c(0,0,1,0,0),2,boxKernel()), c(0,0,0,0,1,1,0,0,0,0))
expect_equal(rescale(c(0,0,1,0,0),2,triangleKernel()), c(0,0,0,0.25,0.75,0.75,0.25,0,0,0))
