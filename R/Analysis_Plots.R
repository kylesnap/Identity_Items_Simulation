library(plyr)
library(tidyverse)
library(arrow)
library(rgl)
library(viridis)
library(showtext)
library(akima)

font_add_google("Open Sans")
showtext_auto()

## Summarising across runs ----

balanced <- open_dataset("./output/balanced")
balanced$schema

computed <- balanced |>
  ungroup() |>
  mutate(
c
  ) |>
  mutate(
    across(
      Pr_XM:Pr_XW, 
      ~ round(. / 0.05) * 0.05, 
      .names = "Bin_{.col}"
    )
  )

summarised <- computed |>
  group_by(Pr_X, D_W, D_X, Bin_Pr_XM, Bin_Pr_XW) |>
  summarise(
    across(
      Bias_Int:Bias_XX,
      list(
        M = ~ mean(., na.rm = TRUE), 
        SD = ~ sd(., na.rm = TRUE),
        FAIL = ~ sum(is.na(.))
      ),
      .names = "{.fn}_{.col}"
    ),
    "Precis_M" = mean(MM / (MM + XM)),
    "Precis_W" = mean(WW / (WW + XW)),
    "Recall_X" = mean(XX / (XX + XM + XW)),
    "Count" = n(),
    .groups = "drop"
  ) |>
  collect() |>
  mutate(
    Pr_X = factor(Pr_X, levels = c(0.1, 0.2, 0.3)),
    D_W = factor(D_W, levels = c(0, 0.5, 0.8)),
    D_X = factor(D_X, levels = c(0, 0.5, 0.8))
  )

## Gridding data to plots ----

make_grid <- function(x, y, z, name, n = 100) {
  #domain <- seq(0.50, 1, by = interval)
  #range <- seq(0.50, 1, by = interval)
  res <- interp(x, y, z, nx = n, ny = n) |>
    interp2xyz(data.frame = TRUE)
  names(res) <- c("Precis_M", "Precis_W", name)
  res
}

gridded <- summarised |>
  ungroup() |>
  group_by(Pr_X, D_W, D_X) |>
  select(Precis_M, Precis_W, M_Bias_Z, M_Bias_XW, M_Bias_XX) |>
  nest(
    "Raw_XZ" = c(Precis_M, Precis_W, M_Bias_Z),
    "Raw_XW" = c(Precis_M, Precis_W, M_Bias_XW),
    "Raw_XX" = c(Precis_M, Precis_W, M_Bias_XX)
  ) |>
  mutate(
    "Interp_XZ" = map(
      Raw_XZ, 
      \(x) make_grid(x$Precis_M, x$Precis_W, x$M_Bias_Z, "Bias_Z")
    ),
    "Interp_XW" = map(
      Raw_XW, 
      \(x) make_grid(x$Precis_M, x$Precis_W, x$M_Bias_XW, "Bias_XW")
    ),
    "Interp_XX" = map(
      Raw_XX, 
      \(x) make_grid(x$Precis_M, x$Precis_W, x$M_Bias_XX, "Bias_XX")
    ),
    .keep = "unused",
  ) |>
  mutate(
    "Interps" = pmap(
      list(Interp_XZ, Interp_XW, Interp_XX),
      \(x, y, z) {
        full_join(
          x, 
          (full_join(y, z, by = c("Precis_M", "Precis_W"))), 
          by = c("Precis_M", "Precis_W")
        ) |>
        drop_na()
      }
    ),
    .keep = "unused"
  ) |>
  unnest(Interps)

count_df <- summarised |>
  group_by(Pr_X, D_W, D_X) |>
  mutate(
    Precis_M,
    Precis_W,
    Prop_XX = case_when(
      FAIL_Bias_XX == 0 ~ "Success",
      (FAIL_Bias_XX / Count) < 0.05 ~ "Partial",
      (FAIL_Bias_XX / Count) > 0.05 ~ "Fail"
    ),
    Prop_XW = case_when(
      FAIL_Bias_XW == 0 ~ "Success",
      (FAIL_Bias_XW / Count) < 0.05 ~ "Partial",
      (FAIL_Bias_XW / Count) > 0.05 ~ "Fail"
    ),
    .keep = "none"
  )

### Plotting: Bias in beta_x, beta_w, beta_z, and beta_z ----

base_df <- gridded |>
  filter(D_W == 0, D_X == 0, Pr_X == 0.3)

base_count <- count_df |>
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

ggplot(base_df, aes(x = Precis_M, y = Precis_W)) +
  geom_tile(aes(fill = Bias_XX, width = 0.01, height = 0.01)) +
  geom_point(data = base_count, size = 2, aes(x = Precis_M, y = Precis_W, shape = Prop_XX)) +
  bias_palette(expression(paste('Bias ', beta["x|m"]))) +
  scale_x_reverse("Precision (M)") +
  scale_y_continuous("Precision (W)") +
  scale_shape_manual(
    "Cell Outcome",
    breaks = c("Success", "Partial", "Fail"),
    labels = c("Success", "<5% Failures", ">5% Failures"),
    values = c(4, 13, 1)
  ) +
  labs(
    title = "Bias in X v. M regression parameter across item precision",
    caption = "Prop. X = 0.3; Cov(X, Z) = Cov(W, Z) = 0"
  ) +
  theme_linedraw(base_family = "Open Sans", base_size = 18) -> bias_xx_plot

ggplot(base_df, aes(x = Precis_M, y = Precis_W)) +
  geom_tile(aes(fill = Bias_XW, width = 0.01, height = 0.01)) +
  geom_point(data = base_count, size = 2, aes(x = Precis_M, y = Precis_W, shape = Prop_XW)) +
  bias_palette(expression(paste('Bias ', beta["w|m"]))) +
  scale_x_reverse("Precision (M)") +
  scale_y_continuous("Precision (W)") +
  scale_shape_manual(
    "Cell Outcome",
    breaks = c("Success", "Partial", "Fail"),
    labels = c("Success", "<5% Failures", ">5% Failures"),
    values = c(4, 13, 1)
  ) +
  labs(
    title = "Bias in W v. M regression parameter across item precision",
    caption = "Prop. X = 0.3; Cov(X, Z) = Cov(W, Z) = 0",
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

delta_df <- gridded |>
  filter(Pr_X == 0.3)

ggplot(delta_df, aes(x = Precis_M, y = Precis_W)) +
  geom_tile(aes(fill = Bias_Z, width = 0.01, height = 0.01)) +
  bias_palette(expression(paste('Bias ', beta["z"]))) +
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
    title = "Bias in Z regression parameter across item precision",
    caption = "Prop. X = 0.3"
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

## Summary of bias over time ----

by_recall <- computed |>
  filter(D_W == 0, D_X == 0) |>
  mutate(
    Pr_X,
    "Recall_X" = XX / (XX + XM + XW),
    "Bias_XX" = OBS_XX - ACT_XX,
    .keep = "none"
  ) |>
  collect() |>
  mutate(Pr_X = factor(Pr_X, levels = c(0.1, 0.2, 0.3)))

alt_background <- theme_linedraw(base_family = "Open Sans", base_size = 18)
alt_background$plot.background <- element_rect(
  fill = alpha("#97d3e9", 0.25), 
  colour = alpha("#97d3e9", 0.25)
)
alt_background$legend.background <- element_rect(fill = alpha("#97d3e9", 0.25))

ggplot(by_recall, aes(Recall_X, Bias_XX, colour = Pr_X)) +
  stat_summary_bin(breaks = seq(1, 0, -0.05)) +
  scale_x_continuous("Recall (X)") +
  scale_y_continuous(expression(paste('Bias ', beta["x|m"]))) +
  scale_colour_brewer("Pr(X)", palette = "Set2") +
  labs(
    title = "Bias in X v. M regression parameter over response recall",
    caption = "Cov(X, Z) = Cov(W, Z) = 0"
  ) +
  alt_background -> plot_recall

ggsave(
  "./output/figure_4_recall.eps",
  plot_recall,
  width = 9,
  height = 9,
  units = "in"
)
