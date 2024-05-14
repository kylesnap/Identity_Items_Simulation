library(profvis)
source("./R/simulation_functions.R")
# Test

cell <- make_cell(
    n = 10,
    probs = c(0.25, 0.25),
    beta = c(0, 0, 1, 0),
    delta = c(0, 0),
    sigma = 1
  )

pmat <- make_pmat(ep_M = c(0, 0), ep_W = c(0.5, 0.5), ep_X = c(0, 0))

result <- run_cell(cell, pmat, reps = 2000)

# Sim table