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
      sample_size = 1000, 
      probabilities = c(
        "M" = (1 - prx)/2, 
        "W" = (1 - prx)/2, 
        "X" = prx
      ),
      true_beta = c(Int = 0, Z = 1, XW = 1, XX = -1),
      deltas = c(W = dw, X = dx),
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
        row_m = c(1, 0, 0), 
        row_w = c(0, 1, 0),
        row_x = c(xm, xw, 1 - (xm + xw))
      )
    ),
    .keep = "none"
  )

three_sim_runs <- cross_join(balanced_cells, balanced_pmats) |>
  mutate(
    N = map_int(Cell, "sample_size"),
    Pr_M = map_dbl(Cell, "Pr_M"),
    Pr_W = map_dbl(Cell, "Pr_W"),
    Pr_X = map_dbl(Cell, "Pr_X"),
    Beta_Int = map_dbl(Cell, "B_Int"),
    Beta_Z = map_dbl(Cell, "B_Z"),
    Beta_XW = map_dbl(Cell, "B_XW"),
    Beta_XX = map_dbl(Cell, "B_XX"),
    D_M = map_dbl(Cell, "D_M"),
    D_W = map_dbl(Cell, "D_W"),
    D_X = map_dbl(Cell, "D_X"),
    Sigma = map_dbl(Cell, "sigma")
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
        run_cell(cell = x, p_mat = y, reps = 1000)
      }
    ),
    .keep = "unused",
    .after = everything()
    ) |>
    unnest(Result)
  result |>
    group_by(Pr_X) |>
    write_dataset(path = "./output/balanced")
}, enable = TRUE, delay_stdout = TRUE, delay_conditions = "condition")
