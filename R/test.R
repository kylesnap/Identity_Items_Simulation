library(profvis)
source("./R/simulation_functions.R")
# Test

cell1 <- make_cell(
  sample_size = 100,
  probabilities = c("M" = 0.5, "W" = 0.3, "X" = 0.2)
)

row_M <- c(0.2, 0.5, 0.3)
row_W <- c(0.1, 0.6, 0.3)
row_X <- c(0.3, 0.4, 0.3)
pmat1 <- make_pmat(row_M, row_W, row_X)

result <- run_cell(cell = cell1, p_mat = pmat1, reps = 10)
