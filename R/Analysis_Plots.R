library(plyr)
library(tidyverse)
library(arrow)
library(showtext)
library(akima)

font_add_google("Open Sans")
showtext_auto()

## Summarising across runs ----

balanced <- open_dataset("./output/balanced")
balanced$schema

computed <- balanced |>
  mutate(
    Pr_X, D_W, D_X,
    Pr_XM = XM / (XM + XW + XX),
    Pr_XW = XW / (XM + XW + XX),
    Bias_Int = OBS_INT - ACT_INT,
    Bias_Z = OBS_Z - ACT_Z,
    Bias_XW = OBS_XW - ACT_XW,
    Bias_XX = OBS_XX - ACT_XX,
  ) |>
  mutate(
    across(
      Pr_XM:Pr_XW, 
      ~ round(. / 0.05) * 0.05, 
      .names = "Bin_{.col}"
    )
  ) |>
  collect()

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

make_grid <- function(x, y, z, name, n = 500) {
  res <- interp(
    x, y, z, 
    nx = n, ny = n,
    duplicate = "mean"
  ) |>
    interp2xyz(data.frame = TRUE)
  names(res) <- c("Precis_M", "Precis_W", name)
  res
}

gridded <- summarised |>
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
    breaks = seq(-1, 1, 0.1)
  )
}

ggplot(base_df, aes(x = Precis_M, y = Precis_W, fill = Bias_XX)) +
  geom_tile() +
  bias_palette(expression(paste('Bias ', beta["x|m"]))) +
  scale_x_continuous("Precision (M)", limits = c(0.55, 1)) +
  scale_y_continuous("Precision (W)", limits = c(0.55, 1)) +
  annotate(
    "label",
    x = 0.55,
    y = 0.55,
    hjust = "left",
    fill = "#dadde8",
    label = "Prop. X = 0.3; Cov(X, Z) = Cov(W, Z) = 0"
  ) -> bias_xx_plot

ggplot(base_df, aes(x = Precis_M, y = Precis_W, fill = Bias_XW)) +
  geom_tile() +
  bias_palette(expression(paste('Bias ', beta["x|m"]))) +
  scale_x_continuous("Precision (M)", limits = c(0.55, 1)) +
  scale_y_continuous("Precision (W)", limits = c(0.55, 1)) +
  annotate(
    "label",
    x = 0.55,
    y = 0.55,
    hjust = "left",
    fill = "#dadde8",
    label = "Prop. X = 0.3; Cov(X, Z) = Cov(W, Z) = 0"
  ) -> bias_xw_plot

delta_df <- gridded |>
  filter(Pr_X == 0.3)

ggplot(delta_df, aes(x = Precis_M, y = Precis_W)) +
  geom_tile(aes(fill = Bias_Z)) +
  bias_palette(expression(paste('Bias ', beta["z"]))) +
  scale_x_reverse("Precision (M)") +
  scale_y_continuous("Precision (W)") +
  facet_grid(
    rows = vars(D_W),
    cols = vars(D_X),
    labeller = labeller(
      .rows = \(x) paste0("Cov(W,Z)=", x),
      .cols = \(x) paste0("Cov(X,Z)=", x)
    )
  ) #-> bias_z_plot

## Summary of bias over time ----

by_recall <- computed |>
  filter(D_W == 0, D_X == 0) |>
  mutate(
    Pr_X,
    "Recall_X" = XX / (XX + XM + XW),
    "SQ_E_XX" = (OBS_XX - ACT_XX)^2,
    .keep = "none"
  ) |>
  collect() |>
  mutate(Pr_X = factor(Pr_X, levels = c(0.1, 0.2, 0.3)))

ggplot(drop_na(by_recall, SQ_E_XX), aes(Recall_X, SQ_E_XX, colour = Pr_X)) +
  stat_summary_bin(
    breaks = seq(1, 0, length.out = 10),
    geom = "point", 
    fun = mean
  ) +
  stat_summary_bin(
    breaks = seq(1, 0, length.out = 10),
    geom = "line", 
    fun = mean
  ) +
  scale_x_reverse("Recall (X)", breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(expression(paste('Mean Sq. Error ', beta["x|m"]))) +
  scale_colour_brewer(
    NULL,
    labels = \(x) paste(expression(Pr(X)), "=", x),
    palette = "Accent"
  ) +
  annotate(
    "label",
    x = 0.01,
    y = 0.01,
    hjust = "right",
    fill = "#dadde8",
    label = "Cov(X,Z) = Cov(W,Z) = 0"
  ) -> plot_recall

## Save images ----

basic_theme <- theme_classic(base_family = "Open Sans", base_size = 18)
#caption_pos <- theme(plot.caption.position = "panel")
dark_background <- basic_theme + theme(
  plot.background = element_rect(fill = "#00053d", colour = NA),
  panel.background = element_rect(fill = "#00053d", colour = NA),
  legend.background = element_rect(fill = "#00053d", colour = NA),
  legend.title = element_text(colour = "#FFFFFF"),
  legend.text = element_text(colour = "#FFFFFF"),
  axis.title = element_text(colour = "#FFFFFF"),
  axis.text = element_text(colour = "#dadde8"),
  axis.line = element_line(colour = "#dadde8"),
  panel.grid.major = element_line(colour = "#a5a9b6"),
  panel.grid.minor = element_line(colour = "#737888"),
)
plot_recall + dark_background

ggsave(
  "./figure/png/figure_1.png", 
  plot_recall + dark_background,
  width = 14.95, 
  height = 8.8, 
  units = "in"
)

ggsave(
  "./figure/eps/figure_1.eps", 
  plot_recall + dark_background,
  width = 14.95, 
  height = 8.8, 
  units = "in"
)

light_background <- basic_theme + theme(
  plot.background = element_rect(fill = "#97d3e9", colour = NA),
  panel.background = element_rect(fill = "#EBEBEA", colour = NA),
  legend.background = element_rect(fill = "#97d3e9", colour = NA),
  panel.grid.major = element_line(colour = "#FFFFFE"),
  strip.background = element_rect(
    fill = "#84C0D6",
    linewidth = 0
  ),
  aspect.ratio = 1
)

ggsave(
  "./figure/png/figure_2.png", 
  bias_xx_plot + light_background,
  width = 15, 
  height = 6.4, 
  units = "in"
)

ggsave(
  "./figure/eps/figure_2.eps", 
  bias_xx_plot + light_background,
  width = 15, 
  height = 6.4, 
  units = "in"
)

bias_xw_plot + light_background

ggsave(
  "./figure/png/figure_3.png", 
  bias_xw_plot + light_background,
  width = 15, 
  height = 6.4, 
  units = "in"
)

ggsave(
  "./figure/eps/figure_3.eps", 
  bias_xw_plot + light_background,
  width = 15, 
  height = 6.4, 
  units = "in"
)

ggsave(
  "./figure/png/figure_4.png", 
  bias_z_plot + light_background,
  width = 15, 
  height = 6.4, 
  units = "in"
)

ggsave(
  "./figure/eps/figure_4.eps", 
  bias_z_plot + light_background,
  width = 15, 
  height = 6.4,
  units = "in"
)