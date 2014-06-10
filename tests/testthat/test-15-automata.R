source(system.file("tests", "testthat", "matches_file.R", package="mmand"))

context("Running cellular automaton simulations")

test_that("cellular automaton simulation works", {
    expect_that(gosperGliderGun(), matches_file("glider_gun_t0.rds"))
    expect_that(gameOfLife(init=gosperGliderGun(),steps=1), matches_file("glider_gun_t1.rds"))
    expect_that(gameOfLife(init=gosperGliderGun(),steps=5), matches_file("glider_gun_t5.rds"))
})
