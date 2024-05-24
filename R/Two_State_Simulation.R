# Two Groups
library(tidyverse)
library(progressr)
library(profvis)
library(beepr)
library(arrow)
source("./R/simulation_functions.R")

two_cells <- expand_grid(
  "N" = c(50, 200, 400, 1000), 
  "P_F" = c(0.5, 0.25, 0.1),
  "Delta_F" = c(0, 0.5, 0.8), 
  "Sigma" = c(0.5, 1, 1.5)
) |>
  mutate(
    "Cell" = pmap(
      list(N, P_F, Delta_F, Sigma),
      \(n, p_f, d_f, sigma) 
      make_cell(n, c(p_f, 0), c(0, 1, 1, 0), c(d_f, 0), sigma)
    )
  )

two_pmats <- expand_grid(
  "Pr_MW" = seq(0, 0.9, 0.05),
  "Pr_WM" = seq(0, 0.9, 0.05)
) |>
  mutate(
    "PMat" = pmap(
      list(Pr_MW, Pr_WM),
      \(mw, wm) make_pmat(c(mw, 0), c(wm, 0), c(0, 0))
    )
  )

two_sim_runs <- cross_join(two_cells, two_pmats)

handlers("debug")
with_progress({
  cl <- makeForkCluster(nnodes = 10)
  registerDoParallel(cl)
  p <- progressor(nrow(two_sim_runs))
  result <- two_sim_runs |>
    mutate("Result" = map2(
      Cell,
      PMat,
      \(x, y) {
        p()
        run_cell(x, y, 500)
      }
    ),
    .keep = "unused"
    ) |>
    unnest(Result)
  write_parquet(result, "./output/raw_two_cells.parquet")
}, enable = TRUE, delay_stdout = TRUE, delay_conditions = "condition")
