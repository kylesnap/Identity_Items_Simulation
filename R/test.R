# Test

cell <- make_cell(
    n = 100,
    probs = c(0.5, 0),
    beta = c(0, 0, 1, 0),
    delta = c(0, 0),
    sigma = 1
  )

pmat <- make_pmat(ep_M = c(0, 0), ep_W = c(0, 0), ep_X = c(0, 0))

print(make_data(cell, pmat))

run_cell(cell, pmat, reps = 1)

# Sim table