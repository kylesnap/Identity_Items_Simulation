library(profvis)
source("./R/simulation_functions.R")
# Test

cell1 <- make_cell(1000, 0.2, 0.3)

row_M <- c(0.2, 0.5, 0.3)
row_W <- c(0.1, 0.6, 0.3)
row_X <- c(0.3, 0.4, 0.3)
pmat1 <- make_pmat(row_M, row_W, row_X)

result <- run_cell(cell1, pmat1, reps = 1, verbose = TRUE)