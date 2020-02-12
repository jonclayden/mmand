# Cellular automaton simulation
expect_equal(gosperGliderGun(), readRDS("glider_gun_t0.rds"))
expect_equal(gameOfLife(init=gosperGliderGun(),steps=1), readRDS("glider_gun_t1.rds"))
expect_equal(gameOfLife(init=gosperGliderGun(),steps=5), readRDS("glider_gun_t5.rds"))
expect_silent(gameOfLife(init=gosperGliderGun(),size=c(40,40),steps=1,viz=TRUE,tick=0))
