library(tidyverse)
library(checkmate)
library(rlang)
library(foreach)
library(doRNG)
library(doParallel)

#' Make a simulation scenario
#' 
#' @param sample_size The number of observations to simulate
#' @param probabilities A named numeric vector of probabilities for each category
#' @param true_beta A named numeric vector of true coefficients
#' @param deltas A named numeric vector of category means
#' @param sigma The standard deviation of the error term
#' 
#' @details The inputs must be named as follows:
#' - `probabilities` must be named `Pr_M`, `Pr_W`, and `Pr_X`
#' - `true_beta` must be named `B_Int`, `B_Z`, `B_XW`, and `B_XX`
#' - `deltas` must be named `D_W` and `D_X` 
#' 
#' @return A list of simulation parameters
make_cell <- function(sample_size,
                      probabilities,
                      true_beta = c("Int" = 0, "Z" = 0, "XW" = 0, "XX" = 0),
                      deltas = c("W" = 0, "X" = 0),
                      sigma = 1) {
  assert_numeric(
    probabilities,
    lower = 0,
    upper = 1,
    len = 3,
    any.missing = FALSE,
    names = "named"
  )
  assert(sum(probabilities) == 1)
  probabilities <- set_names(probabilities, \(x) paste0("Pr_", x))

  assert_numeric(
    true_beta,
    len = 4,
    any.missing = FALSE,
    names = "named"
  )
  true_beta <- set_names(true_beta, \(x) paste0("B_", x))

  assert_numeric(
    deltas,
    len = 2,
    any.missing = FALSE,
    names = "named"
  )
  deltas <- set_names(deltas, \(x) paste0("D_", x))

  assert_count(sample_size)
  assert_number(sigma, lower = .Machine$double.xmin)

  dm <- -((deltas["D_W"] * probabilities["Pr_W"] + 
          deltas["D_X"] * probabilities["Pr_X"]) / 
          (1 - probabilities["Pr_W"] - probabilities["Pr_X"]))[[1]]
  list2(
    "sample_size" = assert_count(sample_size),
    !!!probabilities,
    !!!true_beta,
    "D_M" = dm,
    !!!deltas,
    "sigma" = sigma
  )
}

#' Make a transition matrix
#' 
#' @param row_m,row_w,row_x A numeric vector of probabilities for each category
#' 
#' @return A 3x3 matrix of probabilities
make_pmat <- function(row_m, row_w, row_x) {
  assert_numeric(
    row_m,
    lower = 0,
    upper = 1,
    len = 3,
    any.missing = FALSE
  )
  assert_numeric(
    row_w,
    lower = 0,
    upper = 1,
    len = 3,
    any.missing = FALSE
  )
  assert_numeric(
    row_x,
    lower = 0,
    upper = 1,
    len = 3,
    any.missing = FALSE
  )
  assert(sum(row_m) == 1)
  assert(sum(row_w) == 1)
  assert(sum(row_x) == 1)
  matrix(
    c(row_m, row_w, row_x),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("X=M", "X=W", "X=X"), c("V=M", "V=W", "V=X"))
  )
}

#' Make a sample of continuous variable observations with categorical covariance
#' 
#' @param categories A numeric vector of category labels
#' @param deltas A numeric vector of category means
#' 
#' @details The outcoming variable will always have a standard deviation of 1,
#' and a mean of zero (provided it is called from a cell processed by make_cell
#' 
#' @return A numeric vector of observations of the covariate
make_z <- function(categories, deltas) {
  Z <- numeric(length(categories))
  Z[categories == 1] <- rnorm(
    n = sum(categories == 1), 
    mean = deltas[1], 
    sd = 1
  )
  Z[categories == 2] <- rnorm(
    n = sum(categories == 2),
    mean = deltas[2],
    sd = 1
  )
  Z[categories == 3] <- rnorm(
    n = sum(categories == 3), 
    mean = deltas[3],
    sd = 1
  )
  return(Z)
}

#' Make a sample of observations from a categorical variable and its indicator variable
#' 
#' @param sample_size The number of observations to simulate
#' @param probabilities A numeric vector of probabilities for each category
#' @param p_mat A 3x3 matrix of transition probabilities
#' 
#' @return A list:
#' - X: A numeric vector of category labels
#' - V: A numeric vector of indicator labels
#' - counts: A named numeric vector of counts for each category and indicator combination
make_xv <- function(sample_size, probabilities, p_mat) {
  X <- sample.int(3, size = sample_size, replace = TRUE, prob = probabilities)
  V <- integer(sample_size)
  V[X == 1] <- sample.int(
    3, 
    size = sum(X == 1), 
    replace = TRUE, 
    prob = p_mat[1, ]
  )
  V[X == 2] <- sample.int(
    3, 
    size = sum(X == 2), 
    replace = TRUE, 
    prob = p_mat[2, ]
  )
  V[X == 3] <- sample.int(
    3, 
    size = sum(X == 3), 
    replace = TRUE, 
    prob = p_mat[3, ]
  )

  index <- (X - 1) * 3 + V
  tabs <- tabulate(index, nbins = 9)
  names(tabs) <- c(
    "MM", "MW", "MX",
    "WM", "WW", "WX",
    "XM", "XW", "XX"
  )
  list("X" = X, "V" = V, "counts" = tabs)
}

#' Run a single repetition of a cell with given transition probabilities
#' 
#' @param cell A list of simulation parameters
#' @param p_mat A 3x3 matrix of transition probabilities
#' @param design_mat A matrix of design variables
#' @param Y A numeric vector of responses
#' 
#' @return A list of actual and observed model coefficients and category counts
run_rep <- function(cell, p_mat, design_mat, Y) {
  categories <- make_xv(cell$sample_size, c(cell$Pr_M, cell$Pr_W, cell$Pr_X), p_mat)
  colnames(design_mat) <- c("ACT_Int", "ACT_Z", "ACT_XW", "ACT_XX")
  design_mat[, 2] <- make_z(categories$X, c(cell$D_M, cell$D_W, cell$D_X))
  design_mat[, 3] <- categories$X == 2
  design_mat[, 4] <- categories$X == 3
  Y <- (design_mat %*% c(cell$B_Int, cell$B_Z, cell$B_XW, cell$B_XX)) + 
    rnorm(cell$sample_size, sd = cell$sigma)
  mod_act <- lm.fit(design_mat, Y)$coefficients
  colnames(design_mat) <- c("OBS_INT", "OBS_Z", "OBS_XW", "OBS_XX")
  design_mat[, 3] <- categories$V == 2
  design_mat[, 4] <- categories$V == 3
  mod_obs <- lm.fit(design_mat, Y)$coefficients
  list2(
    !!!mod_act,
    !!!mod_obs,
    !!!categories$counts
  )
}

#' Run a set number of repetitions of a cell with given transition probabilities
#' 
#' @param cell A list of simulation parameters
#' @param p_mat A 3x3 matrix of transition probabilities
#' @param reps The number of repetitions to run
#' 
#' @details The object produced by this function can get incredibly long; consider
#' pre-processing repetition results here or store the results in a binary format
#' 
#' @return A tibble of results, containing the actual and observed model 
#' coefficients and category counts
run_cell <- function(cell, p_mat, reps = 1) {
  print("cell")
  Y <- rep(0, cell$sample_size)
  design_mat <- matrix(
    c(rep(1, cell$sample_size), rep(0, cell$sample_size * 3)), 
    nrow = cell$sample_size
  )
  loop <- foreach(
    1:reps,
    .options.RNG = 420
  )
  result <- if (getDoParRegistered()) {
    loop %dorng% run_rep(cell, p_mat, design_mat, Y)
  } else {
    loop %do% run_rep(cell, p_mat, design_mat, Y)
  }
  map(result, \(x) as_tibble_row(x)) |>
    list_rbind()
}
 