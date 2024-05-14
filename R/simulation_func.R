library(tidyverse)
library(checkmate)
library(rlang)
library(cpp11)
cpp_source("./src/simulation.cpp")

#' Create and validate a cell object for simulation
#'
#' @param n The number of observations.
#' @param probs A numeric vector of length 2 containing the probabilities for each category.
#' @param beta A numeric vector of length 4 representing beta coefficients.
#' @param delta A numeric vector of length 2 representing delta coefficients.
#' @param sigma The standard deviation.
#' @return A list representing the cell object.
make_cell <- function(n, probs, beta, delta, sigma) {
  assert_numeric(
    probs,
    lower = 0,
    upper = 1,
    len = 2,
    any.missing = FALSE
  )
  assert_numeric(beta, len = 4)
  assert_numeric(delta, len = 2)
  assert(probs[1] + probs[2] <= 1)
  list(
    "n" = assert_count(n),
    "probs" = c(
      "Pr_M" = 1 - probs[1] - probs[2],
      "Pr_W" = probs[1],
      "Pr_X" = probs[2]
    ),
    "beta" = c(
      "int" = beta[1],
      "Z" = beta[2],
      "XW" = beta[3],
      "XX" = beta[4]
    ),
    "delta" = c(
      "D_M" = -(delta[1] * probs[1] + delta[2] * probs[2]) / (1 - probs[1] - probs[2]),
      "D_W" = delta[1],
      "D_X" = delta[2]
    ),
    "sigma" = assert_number(sigma, lower = .Machine$double.xmin)
  )
}

#' Create a transition probability matrix
#'
#' @param ep_M The probabilities of transitioning from M.
#' @param ep_W The probabilities of transitioning from W.
#' @param ep_X The probabilities of transitioning from X.
#' @return A transition probability matrix.
make_pmat <- function(ep_M, ep_W, ep_X) {
  assert_numeric(
    ep_M,
    lower = 0,
    upper = 1,
    len = 2,
    any.missing = FALSE
  )
  assert_numeric(
    ep_W,
    lower = 0,
    upper = 1,
    len = 2,
    any.missing = FALSE
  )
  assert_numeric(
    ep_X,
    lower = 0,
    upper = 1,
    len = 2,
    any.missing = FALSE
  )
  assert(sum(ep_M) <= 1)
  assert(sum(ep_W) <= 1)
  assert(sum(ep_X) <= 1)
  matrix(
    c(1 - sum(ep_M),  ep_M[1], ep_M[2],
      ep_W[1], 1 - sum(ep_W), ep_W[2],
      ep_X[1], ep_X[2], 1 - sum(ep_X)
    ),
    nrow = 3,
    byrow = TRUE
  )
}

#' Invent Z and Betas values
make_z <- function(X, deltas) {
  from_M <- rnorm(length(X), mean = deltas["D_M"], sd = 1)
  from_W <- rnorm(length(X), mean = deltas["D_W"], sd = 1)
  from_X <- rnorm(length(X), mean = deltas["D_X"], sd = 1)
  from_M * (X == 1) + from_W * (X == 2) + from_X * (X == 3)
}

#' Mis-categorize observations
#'
#' @param X The original categories.
#' @param pmat The transition probability matrix.
#' @return A vector of mis-categorized observations.
miscategorize <- function(X, pmat) {
  from_M <- sample.int(3, size = length(X), replace = TRUE, prob = pmat[1,])
  from_W <- sample.int(3, size = length(X), replace = TRUE, prob = pmat[2,])
  from_X <- sample.int(3, size = length(X), replace = TRUE, prob = pmat[3,])
  from_M * (X == 1) + from_W * (X == 2) + from_X * (X == 3)
}

#' Create data based on cell and transition probability matrix
#'
#' @param cell The cell object.
#' @param pmat The transition probability matrix.
#' @return A list containing the actual and observed matrices, Y values, and counts.
make_data <- function(cell, pmat) {
  XV <- make_x_v_cpp(cell$n, cell$probs, pmat[1,], pmat[2,], pmat[3,])
  Z <- make_z(XV$X, cell$delta)
  mat_act <- matrix(c(rep(1, cell$n), Z, XV$X == 2, XV$X == 3), nrow = cell$n)
  mat_obs <- matrix(c(rep(1, cell$n), Z, XV$V == 2, XV$V == 3), nrow = cell$n)
  colnames(mat_act) <- c("INT", "Z", "X_W", "X_X")
  colnames(mat_obs) <- c("INT", "Z", "V_W", "V_X")
  Y <- (mat_act %*% cell$beta) + rnorm(cell$n, sd = sqrt(cell$sigma))
  models <- fit_models_cpp(mat_act, mat_obs, Y);
  count <- count_cases_cpp(XV$X, XV$V)
  # X <- sample.int(3, size = cell$n, replace=TRUE, prob = cell$probs)
  # Z <- make_z(X, cell$delta)
  # mat_act <- matrix(c(rep(1, cell$n), Z, X == 2, X == 3), nrow = cell$n)
  # colnames(mat_act) <- c("INT", "Z", "X_W", "X_X")
  # V <- miscategorize(X, pmat)
  # mat_obs <- matrix(c(rep(1, cell$n), Z, V == 2, V == 3), nrow = cell$n)
  # colnames(mat_obs) <- c("INT", "Z", "V_W", "V_X")
  # Y <- (mat_act %*% cell$beta) + rnorm(cell$n, sd = sqrt(cell$sigma))
  # list(
  #   "actual" = mat_act,
  #   "observed" = mat_obs,
  #   "Y" = Y,
  #   count = list(
  #     "X" = X,
  #     "V" = V
  #   )
  # )
}

#' Run linear regression models
#'
#' @param data A list containing actual and observed data matrices and response variable.
#' @return A list containing coefficients and residual sum of squares for actual and observed models.
run_models <- function(data) {
  print(class(data$actual))
  print(class(data$observed))
  act <- fastLmPure(
    X = data$actual,
    y = data$Y
  )
  obs <- fastLmPure(
    X = data$observed,
    y = data$Y
  )
  count <- c(
    "MM" = sum(data$count$X == 1 & data$count$V == 1),
    "MW" = sum(data$count$X == 1 & data$count$V == 2),
    "MX" = sum(data$count$X == 1 & data$count$V == 3),
    "WM" = sum(data$count$X == 2 & data$count$V == 1),
    "WW" = sum(data$count$X == 2 & data$count$V == 2),
    "WX" = sum(data$count$X == 2 & data$count$V == 3),
    "XM" = sum(data$count$X == 3 & data$count$V == 1),
    "XW" = sum(data$count$X == 3 & data$count$V == 2),
    "XX" = sum(data$count$X == 3 & data$count$V == 3)
  )
  names(act$coefficients) <- paste("ACT", names(act$coefficients), sep = "_")
  names(obs$coefficients) <- paste("OBS", names(obs$coefficients), sep = "_")
  list2(
    !!!act$coefficients,
    !!!obs$coefficients,
    !!!count,
    "ACT_rank" = act$rank,
    "OBS_rank" = obs$rank
    #"actual_RSS" = sum(act$residuals^2),
    #"observed_RSS" = sum(obs$residuals^2)
  )
}

#' Run simulations for a single cell
#'
#' @param cell The cell object.
#' @param pmat The transition probability matrix.
#' @param reps The number of repetitions for the cell
#' @return A list of results from the runs
run_cell <- function(cell, pmat, reps = 1) {
  map(
    c(1:reps),
    \(i) run_models(make_data(cell, pmat))
  )
}