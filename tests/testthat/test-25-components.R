context("Connected components")

test_that("the connected components algorithm works", {
    data <- c(0,0,1,0,0,0,1,1,1,0,0)
    kernel <- c(1,1,1)
    
    # Any labelling is valid, but the algorithm is deterministic so should be consistent
    expect_that(symmetric(kernel), is_true())
    expect_that(components(data,kernel), equals(c(NA,NA,2,NA,NA,NA,1,1,1,NA,NA)))
    
    # 1 1 0
    # 1 0 1
    # 0 1 1
    # The two components are connected only diagonally, so the kernel matters
    data <- matrix(c(1,1,0,1,0,1,0,1,1), 3, 3)
    boxKernel <- shapeKernel(c(3,3), type="box")
    diamondKernel <- shapeKernel(c(3,3), type="diamond")
    expect_that(components(data,boxKernel), equals_reference("fan_components_box.rds"))
    expect_that(components(data,diamondKernel), equals_reference("fan_components_diamond.rds"))
    
    # 0 1 1
    # 0 0 0
    # 0 1 1
    # Make sure there is no inappropriate wrap-around
    data <- matrix(c(0,0,0,1,0,1,1,0,1), 3, 3)
    result <- components(data, shapeKernel(c(3,3)))
    expect_that(result[1,3] == result[3,3], is_false())
})
