context("Cellular automaton simulations")

test_that("cellular automaton simulation works", {
    expect_that(gosperGliderGun(), equals_reference("glider_gun_t0.rds"))
    expect_that(gameOfLife(init=gosperGliderGun(),steps=1), equals_reference("glider_gun_t1.rds"))
    expect_that(gameOfLife(init=gosperGliderGun(),steps=5), equals_reference("glider_gun_t5.rds"))
})
