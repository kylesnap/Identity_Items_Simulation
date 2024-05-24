library(tidyverse)
library(arrow)
library(rgl)
library(viridis)

balanced <- open_dataset("./output/balanced_2024-05-24")

computed <- balanced |>
  mutate(
    AF_INT = OBS_INT / ACT_INT,
    AF_Z = OBS_Z / ACT_Z,
    AF_XW = OBS_XW / ACT_XW,
    AF_XX = OBS_XX / ACT_XX
  ) |>
  select(N, Pr_M, Pr_W, Pr_X, Pr_XM, Pr_XW, Pr_XX, D_M, D_W, D_X, Sigma_E, AF_INT:AF_XX)

rec_prec <- computed |>
  collect()

simple_trials <- rec_prec |>
  filter(N == 1000, Delta_X == 0, Delta_W == 0, Sigma == 1) |>
  group_by(Pr_X, Precision_M, Precision_W) |>
  summarise(
    AF_INT = mean(AF_INT, na.rm = TRUE),
    NA_INT = sum(is.na(AF_INT)),
    AF_XW = mean(AF_XW, na.rm = TRUE),
    NA_XW = sum(is.na(AF_XW)),
    AF_XX = mean(AF_XX, na.rm = TRUE),
    NA_XX = sum(is.na(AF_XX)),
    .groups = "drop"
  ) |>
  collect() |>
  mutate(
    NA_INT = case_when(
      NA_INT == 500 ~ "All",
      NA_INT > 0 ~ "Some",
      NA_INT == 0 ~ "None"
    ),
    NA_XW = case_when(
      NA_XW == 500 ~ "All",
      NA_XW > 0 ~ "Some",
      NA_XW == 0 ~ "None"
    ),
    NA_XX = case_when(
      NA_XX == 500 ~ "All",
      NA_XX > 0 ~ "Some",
      NA_XX == 0 ~ "None"
    )
  )

ggplot(simple_trials, aes(x = Precision_M, y = Precision_W, colour = AF_XX, shape = NA_XX)) +
  geom_point(size = 2) +
  scale_x_reverse("Precision (M)") +
  scale_y_continuous("Precision (W)") +
  scale_color_viridis("Attenuation", option = "viridis", na.value = "#bc3754") +
  scale_fill_manual(values = NA) +
  guides(fill = guide_legend("Attenuation ~ Inf", override.aes=list(colour="#bc3754"))) +
  scale_shape_manual(
    "Run Outcome",
    breaks = c("All", "Some", "None"),
    labels = c("No Results", "Partial Results", "Full Results"),
    values = c(4, 13, 16) 
  ) +
  theme_light() +
  facet_grid(cols = vars(Pr_X), labeller = as_labeller(\(x) paste("Pr(X) =", x)))

ggplot(simple_trials, aes(x = Precision_M, y = Precision_W, colour = AF_XX, shape = NA_XW)) +
  geom_point(size = 2) +
  scale_x_reverse("Precision (M)") +
  scale_y_continuous("Precision (W)") +
  scale_color_viridis("Attenuation", option = "viridis", na.value = "#bc3754") +
  scale_fill_manual(values = NA) +
  guides(fill = guide_legend("Attenuation ~ Inf", override.aes=list(colour="#bc3754"))) +
  scale_shape_manual(
    "Run Outcome",
    breaks = c("All", "Some", "None"),
    labels = c("No Results", "Partial Results", "Full Results"),
    values = c(4, 13, 16) 
  ) +
  theme_light() +
  facet_grid(cols = vars(Pr_X), labeller = as_labeller(\(x) paste("Pr(X) =", x)))