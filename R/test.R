library(profvis)
source("./R/simulation_functions.R")
# Test

cell <- make_cell(
    n = 1000,
    probs = c(0., 0.5),
    beta = c(0, 1, 1, 2),
    delta = c(0, 0),
    sigma = 1
  )

pmat <- make_pmat(ep_M = c(0, 0), ep_W = c(0.0, 0.0), ep_X = c(0, 0))

result <- run_cell(cell, pmat, reps = 20)
#})
result_table <- lapply(result, \(x) as_tibble(x)) |>
  list_rbind()

cl <- makePSOCKcluster(5)
clusterCall(cl, function() { source("./R/simulation_functions.R") })
registerDoParallel(cl)
profvis({
  result <- run_cell(cell, pmat, reps = 2000)
})
stopCluster(cl)