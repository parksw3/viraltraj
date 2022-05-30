source("../R/viraltraj.R")
source("../R/simepi.R")

simulation_base_example <- simepi(seed=101)

save("simulation_base_example", file="simulation_base_example.rda")
