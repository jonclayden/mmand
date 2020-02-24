# Cellular automaton simulation
expect_equal_to_reference(gosperGliderGun(), "glider_gun_t0.rds")
expect_equal_to_reference(gameOfLife(init=gosperGliderGun(),steps=1), "glider_gun_t1.rds")
expect_equal_to_reference(gameOfLife(init=gosperGliderGun(),steps=5), "glider_gun_t5.rds")
expect_silent(gameOfLife(init=gosperGliderGun(),size=c(40,40),steps=1,viz=TRUE,tick=0))
