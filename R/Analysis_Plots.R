library(tidyverse)
library(arrow)
library(rgl)
library(viridis)

balanced <- open_dataset("./output/balanced")
balanced$schema

computed <- balanced |>
  mutate(
    "PRECISION_M" = MM / (MM + WM + XM),
    "PRECISION_W" = WW / (WW + MW + XW),
    "RECALL_X" = XX / (XX + XM + XW),
    "DIF_INT" = OBS_INT - ACT_INT,
    "DIF_Z" = OBS_Z - ACT_Z,
    "DIF_XW" = OBS_XW - ACT_XW,
    "DIF_XX" = OBS_XX - ACT_XX
  )

count_na <- \(x) sum(is.na(x))
balanced_summary <- computed |>
  mutate(
    "R_PRECISION_M" = round(PRECISION_M/0.05) * 0.05,
    "R_PRECISION_W" = round(PRECISION_W/0.05) * 0.05,
  ) |>
  group_by(Pr_X, D_W, D_X, R_PRECISION_M, R_PRECISION_W) |>
  summarise(
    across(
      DIF_INT:DIF_XX, 
      list(mean, sd, ~ sum(is.na(.)) / n()),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  collect() |>
  rename_with(
    \(x) str_replace_all(x, c("1$" = "M", "2$" = "SD", "3$" = "FPROP")),
    .cols = DIF_INT_1:DIF_XX_3
  )
                                
  mutate(
    FAILED_XX = FAILED_XX / TOTAL_MOD,
    "Failure_Prop" = case_when(
      FAILED_XX == 1 ~ "All",
      FAILED_XX > 0.5 ~ "Over_Fifty",
      FAILED_XX > 0 ~ "Under_Fifty",
      FAILED_XX == 0 ~ "None"
    )
  )

ggplot(balanced_summary, aes(x = SPE_M_R, y = SPE_W_R, colour = DIF_XX)) +
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

ggplot(simple_trials, aes(x = Precision_M, y = Precision_W, colour = DIF_XX, shape = NA_XW)) +
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