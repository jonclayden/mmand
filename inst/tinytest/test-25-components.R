# Connected components
data <- c(0,0,1,0,0,0,1,1,1,0,0)
kernel <- c(1,1,1)

# Any labelling is valid, but the algorithm is deterministic so should be consistent
expect_true(symmetric(kernel))
expect_equal(components(data,kernel), c(NA,NA,2,NA,NA,NA,1,1,1,NA,NA))

# 1 1 0
# 1 0 1
# 0 1 1
# The two components are connected only diagonally, so the kernel matters
data <- matrix(c(1,1,0,1,0,1,0,1,1), 3, 3)
boxKernel <- shapeKernel(c(3,3), type="box")
diamondKernel <- shapeKernel(c(3,3), type="diamond")
expect_equal_to_reference(components(data,boxKernel), "fan_components_box.rds")
expect_equal_to_reference(components(data,diamondKernel), "fan_components_diamond.rds")

# 0 1 1
# 0 0 0
# 0 1 1
# Make sure there is no inappropriate wrap-around
data <- matrix(c(0,0,0,1,0,1,1,0,1), 3, 3)
result <- components(data, shapeKernel(c(3,3)))
expect_false(result[1,3] == result[3,3])
