# Test

cell <- make_cell(
    n = 100,
    probs = c(0.25, 0.25),
    beta = c(0, 0, 1, 0),
    delta = c(0, 0),
    sigma = 1
  )

pmat <- make_pmat(ep_M = c(0.5, 0.5), ep_W = c(0, 0), ep_X = c(0, 0))

make_x_v_cpp(cell$n, cell$probs, pmat[1,], pmat[2,], pmat[3,])

make_data(cell, pmat)

run_cell(cell, pmat, reps = 1)

# Sim table