library(tidyverse)
library(arrow)
library(rgl)
library(viridis)
library(showtext)

font_add_google("Open Sans")
showtext_auto()

balanced <- open_dataset("./output/balanced")
balanced$schema

computed <- balanced |>
  group_by(Pr_X, D_W, D_X) |>
  mutate(
    N,
    Pr_M = (MM + MW + MX) / N,
    Pr_W = (WW + WM + WX) / N,
    Pr_X = (XX + XM + XW) / N,
    Precision_M = MM / (MM + XM),
    Precision_W = WW / (WW + XW),
    Recall_X = XX / (XX + XM + XW)
  ) |>
  compute()

round_to_step <- function(var, step = 1) {
  round(var/step) * step
}

regrouped <- computed |>
  mutate(
    across(c(Pr_M, Pr_W, Pr_X), ~ round(./0.1) * 0.1),
    across(Precision_M:Precision_W, ~ round(./0.025) * 0.025)
  )

by_precision <- regrouped |>
  group_by(Pr_X, D_W, D_X, Precision_M, Precision_W) |>
  summarise(
    across(ACT_INT:OBS_XX, ~ mean(., na.rm = TRUE)),
    "MSE_INT" = mean((OBS_INT - ACT_INT)^2),
    "MSE_Z" = mean((OBS_Z - ACT_Z)^2),
    "MSE_XW" = mean((OBS_XW - ACT_XW)^2),
    "MSE_XX" = mean((OBS_XX - ACT_XX)^2, na.rm = TRUE),
    "XX_NA" = sum(is.na(OBS_XX)),
    "Count" = n()
  ) |>
  ungroup() |>
  collect() |>
  mutate(
    "BIAS_INT" = OBS_INT - ACT_INT,
    "BIAS_Z" = OBS_Z - ACT_Z,
    "BIAS_XW" = OBS_XW - ACT_XW,
    "BIAS_XX" = OBS_XX - ACT_XX,
    Pr_X = factor(Pr_X, levels = c(0.1, 0.2, 0.3))
  ) |>
  mutate(
    XX_NA_PROP = case_when(
      XX_NA / Count > 0.05 ~ "Failure",
      XX_NA / Count == 0 ~ "Success",
      XX_NA / Count <= 0.05 ~ "Part_Failure"
    )
  )

base_df <- by_precision |>
  filter(D_W == 0, D_X == 0, Pr_X == 0.3)

bias_palette <- function(name, ...) {
  scale_fill_gradient2(
    name,
    low = "#364B9A",
    mid = "#EAECCC",
    high = "#A50026",
    midpoint = 0,
    na.value = "#FFFFFF",
    ...
  )
}

no_outline <- theme_linedraw(base_family = "Open Sans", base_size = 18)
no_outline$plot.background = element_rect(linewidth = 0)

ggplot(base_df, aes(x = Precision_M, y = Precision_W, fill = BIAS_XX, shape = XX_NA_PROP)) +
  geom_point(size = 3) +
  bias_palette(expression(paste('Bias ', beta[x]))) +
  scale_x_reverse("Precision (M)") +
  scale_y_continuous("Precision (W)") +
  scale_shape_manual(
    "Failure Rate",
    labels = c(">5% Trials", "<5% Trials", "No Trials"),
    breaks = c("Failure", "Part_Failure", "Success"),
    values = c(23, 25, 21)
  ) +
  labs(
    title = expression(
      paste("Bias in ", beta[x], " across levels of survey item precision")
    ),
    caption = "Sample Size = 1000; Prop. X = 0.3; Cov(X, Z) = Cov(W, Z) = 0"
    ) +
  theme_linedraw(base_family = "Open Sans", base_size = 18) -> bias_xx_plot

ggplot(base_df, aes(x = Precision_M, y = Precision_W, fill = BIAS_XW)) +
  geom_point(size = 3, shape = 21) +
  bias_palette(expression(paste('Bias ', beta[w]))) +
  scale_x_reverse("Precision (M)") +
  scale_y_continuous("Precision (W)") +
  labs(
    title = expression(
      paste("Bias in ", beta[w], " across levels of survey item precision")
    ),
    caption = "Sample Size = 1000; Prop. X = 0.3; Cov(X, Z) = Cov(W, Z) = 0"
  ) +
  theme_linedraw(base_family = "Open Sans", base_size = 18) -> bias_xw_plot

ggsave(
  "./output/figure_1_xx.eps",
  bias_xx_plot,
  width = 9,
  height = 9,
  units = "in",
  bg = "transparent"
)

ggsave(
  "./output/figure_2_xw.eps",
  bias_xw_plot,
  width = 9,
  height = 9,
  units = "in",
  bg = "transparent"
)

delta_df <- by_precision |>
  filter(Pr_X == 0.3, D_W != 0.8, D_X != 0.8)

ggplot(delta_df, aes(x = Precision_M, y = Precision_W, fill = BIAS_Z)) +
  geom_point(size = 2, shape = 21) +
  bias_palette(expression(paste('Bias ', beta[z]))) +
  scale_x_reverse("Precision (M)") +
  scale_y_continuous("Precision (W)") +
  facet_grid(
    rows = vars(D_W), 
    cols = vars(D_X),
    labeller = labeller(
      .rows = \(x) paste("Cov(W,Z) = ", x),
      .cols = \(x) paste("Cov(X,Z) = ", x)
    )
  ) +
  labs(
    title = expression(
      paste("Bias in ", beta[z], " across levels of survey item precision")
    ),
    caption = "Sample Size = 1000; Prop. X = 0.3"
  ) +
  theme_linedraw(base_family = "Open Sans", base_size = 18) -> bias_z_plot

ggsave(
  "./output/figure_3_z.eps",
  bias_z_plot,
  width = 9,
  height = 9,
  units = "in",
  bg = "transparent"
)

by_recall <- regrouped |>
  filter(D_X == 0, D_W == 0) |>
  collect() |>
  mutate(
    Recall_X_Bin = cut_width(
      Recall_X, 
      0.05, 
      boundary = 0,
      labels = FALSE
    )
  ) |>
  group_by(Pr_X, Recall_X_Bin) |>
  summarise(
    "Center" = max(Recall_X),
    "BIAS_XX" = mean_cl_normal((OBS_XX - ACT_XX)),
    "XX_NA" = sum(is.na(OBS_XX)),
    "Count" = n(),
    .groups = "drop"
  ) |>
  unpack(BIAS_XX) |>
  ungroup() |>
  mutate(
    "Recall_X_Bin" = round(Center/0.05) * 0.05,
    Pr_X = factor(Pr_X, levels = c(0.1, 0.2, 0.3)),
    XX_NA_PROP = case_when(
      XX_NA / Count > 0.05 ~ "Failure",
      XX_NA / Count == 0 ~ "Success",
      XX_NA / Count <= 0.05 ~ "Part_Failure"
    ),
    .keep = "unused"
  )

alt_background <- theme_linedraw(base_family = "Open Sans", base_size = 18)
alt_background$plot.background <- element_rect(fill = alpha("#97d3e9", 0.25), colour = alpha("#97d3e9", 0.25))
alt_background$legend.background <- element_rect(fill = alpha("#97d3e9", 0.25))

ggplot(by_recall, aes(Recall_X_Bin, y = y, ymin = ymin, ymax = ymax, 
                      colour = Pr_X)) +
  geom_linerange(linewidth = 1) +
  geom_point(size = 2, aes(shape = Pr_X)) +
  scale_x_reverse("Recall (X)") +
  scale_y_continuous(expression(paste('Bias ', beta[x]))) +
  scale_colour_viridis_d("Prop. N.B.") +
  labs(
    title = expression(
      paste("Bias of ", beta[z], " over non-binary gender item recall")
    ),
    caption = "Sample Size = 1000; Cov(X, Z) = Cov(W, Z) = 0"
  ) +
  guides(shape = "none") +
  alt_background -> plot_recall

ggsave(
  "./output/figure_4_recall.eps",
  plot_recall,
  width = 9,
  height = 9,
  units = "in"
)
  
