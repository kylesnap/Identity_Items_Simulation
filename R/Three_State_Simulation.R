# Two Groups
library(tidyverse)
library(progressr)
library(arrow)
source("./R/simulation_functions.R")

balanced_cells <- expand_grid(
  "Pr_X" = c(0.1, 0.2, 0.3),
  "Delta_W" = c(0, 0.5, 0.8),
  "Delta_X" = c(0, 0.5, 0.8),
) |>
  mutate(
    "Cell" = pmap(list(Pr_X, Delta_W, Delta_X), 
         \(prx, dw, dx) make_cell(
      N = 1000, 
      Pr_W = round((1 - prx)/2, 2),
      Pr_X = prx, 
      Beta_Z = 1, 
      Beta_XW = 1,
      Beta_XX = -1,
      Delta_W = dw,
      Delta_X = dx
    )),
    .keep = "none"
  )

balanced_pmats <- expand_grid(
  "Pr_XM" = seq(0, 1, 0.05),
  "Pr_XW" = seq(0, 1, 0.05)
) |>
  filter(Pr_XM + Pr_XW <= 1) |>
  mutate(
    "PMat" = map2(
      Pr_XM,
      Pr_XW,
      \(xm, xw) make_pmat(
        row_M = c(1, 0, 0), 
        row_W = c(0, 1, 0),
        row_X = c(xm, xw, 1 - (xm + xw))
      )
    ),
    .keep = "none"
  )

three_sim_runs <- cross_join(balanced_cells, balanced_pmats) |>
  mutate(
    N = map_int(Cell, "n"),
    Pr_M = map_dbl(Cell, \(x) pluck(x, "probs", "Pr_M")),
    Pr_W = map_dbl(Cell, \(x) pluck(x, "probs", "Pr_W")),
    Pr_X = map_dbl(Cell, \(x) pluck(x, "probs", "Pr_X")),
    Beta_Int = map_dbl(Cell, \(x) pluck(x, "beta", "Int")),
    Beta_Z = map_dbl(Cell, \(x) pluck(x, "beta", "Z")),
    Beta_XW = map_dbl(Cell, \(x) pluck(x, "beta", "XW")),
    Beta_XX = map_dbl(Cell, \(x) pluck(x, "beta", "XX")),
    D_M = map_dbl(Cell, \(x) pluck(x, "delta", "D_M")),
    D_W = map_dbl(Cell, \(x) pluck(x, "delta", "D_W")),
    D_X = map_dbl(Cell, \(x) pluck(x, "delta", "D_X")),
    Sigma_E = map_dbl(Cell, "sigma")
  )

handlers("pbcol")
with_progress({
  cl <- makeForkCluster(nnodes = 10)
  registerDoParallel(cl)
  p <- progressor(nrow(three_sim_runs))
  result <- three_sim_runs |>
    mutate("Result" = map2(
      Cell,
      PMat,
      \(x, y) {
        p()
        run_cell(x, y, 1000)
      }
    ),
    .keep = "unused",
    .after = everything()
    ) |>
    unnest(Result)
  print(sapply(result, function(x) length(names(x))))
  result |>
    group_by(Pr_X) |>
    write_dataset(path = "./output/balanced")
}, enable = TRUE, delay_stdout = TRUE, delay_conditions = "condition")
