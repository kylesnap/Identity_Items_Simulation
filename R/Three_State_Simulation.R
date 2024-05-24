# Two Groups
library(tidyverse)
library(progressr)
library(arrow)
source("./R/simulation_functions.R")

three_cells <- expand_grid(
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

three_pmats <- expand_grid(
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

three_sim_runs <- cross_join(three_cells, three_pmats) |>
  mutate(
    Cell_Params = map(Cell, \(x) cell_as_row(x)),
    PMat_Params = map(PMat, \(x) pmat_as_row(x))
  ) |>
  unnest_wider(Cell_Params:PMat_Params)

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
  write_dataset(
    result,
    path = paste0("./output/balanced_", Sys.Date()),
    partitioning = c("Pr_X", "D_W", "D_X")
  )
}, enable = TRUE, delay_stdout = TRUE, delay_conditions = "condition")
