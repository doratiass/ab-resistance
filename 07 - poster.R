# ----------------------------
# Libraries
# ----------------------------
library(tidyverse)
library(glue)
library(showtext)

# ----------------------------
# 1) Data
# ----------------------------
nnh_df <- tribble(
  ~reference,
  ~matched,
  ~resistance_to,
  ~NNH,
  ~RD_pct,
  ~CI_low_pct,
  ~CI_high_pct,
  # ESBL proxy: 3rd-generation ceph resistance
  "1st CEP",
  "AMC",
  "3rd CEP",
  27.29,
  3.7,
  2.0,
  5.3,
  "1st CEP",
  "2nd CEP",
  "3rd CEP",
  11.73,
  8.5,
  6.8,
  10.2,
  "1st CEP",
  "FQ",
  "3rd CEP",
  18.82,
  5.3,
  3.6,
  7.1,
  "AMC",
  "2nd CEP",
  "3rd CEP",
  23.56,
  4.2,
  2.9,
  5.6,
  "FQ",
  "2nd CEP",
  "3rd CEP",
  46.95,
  2.1,
  0.4,
  3.9,
  # Nitrofurantoin resistance
  "1st CEP",
  "AMC",
  "nitrofurantoin",
  44.05,
  2.3,
  0.1,
  4.4,
  "1st CEP",
  "2nd CEP",
  "nitrofurantoin",
  139.52,
  0.7,
  -1.3,
  2.8,
  "1st CEP",
  "FQ",
  "nitrofurantoin",
  -274.12,
  -0.4,
  -2.5,
  1.8,
  "AMC",
  "2nd CEP",
  "nitrofurantoin",
  -201.23,
  -0.5,
  -2.1,
  1.1,
  "FQ",
  "2nd CEP",
  "nitrofurantoin",
  84.88,
  1.2,
  -0.7,
  3.0
) %>%
  mutate(
    pair = glue("{matched} vs {reference}"),
    group = if_else(
      resistance_to == "3rd CEP",
      "ESBL resistance",
      "Nitrofurantoin resistance"
    ),
    ci_spans_zero = CI_low_pct <= 0 & CI_high_pct >= 0,
    # NNH / NNB label only when CI excludes 0 (i.e., direction is clear)
    NNH_num = if_else(RD_pct != 0, 100 / abs(RD_pct), NA_real_),
    NNH_label = case_when(
      !ci_spans_zero & RD_pct > 0 ~ glue("NNH {round(NNH_num)}"),
      !ci_spans_zero & RD_pct < 0 ~ glue("NNB {round(NNH_num)}"),
      TRUE ~ "No clear effect"
    )
  )

# ----------------------------
# 2) Ordering & aesthetics
# ----------------------------
# Order pairs within each outcome by descending RD (so largest harm at top)
nnh_df <- nnh_df %>%
  group_by(resistance_to) %>%
  arrange(desc(RD_pct), .by_group = TRUE) %>%
  mutate(pair = factor(pair, levels = rev(unique(pair)))) %>%
  ungroup()

# Okabe–Ito palette (color-blind friendly)
pal_esbl <- "#71A632" # warm accent
pal_nitro <- "#742D91" # cool accent
pal_grey <- "#8C8C8C" # muted for CIs crossing 0

# Map color by group, then reduce alpha when CI crosses zero
nnh_df <- nnh_df %>%
  mutate(
    color_grp = if_else(group == "ESBL resistance", pal_esbl, pal_nitro),
    alpha_val = if_else(ci_spans_zero, 0.45, 1.0),
    # Where to place label (right if RD>=0, left if RD<0)
    # label_x = if_else(RD_pct >= 0, CI_high_pct + 0.3, CI_low_pct - 0.3),
    label_x = abs(CI_high_pct) + 0.3,
    label_hj = 0
    # label_hj = if_else(RD_pct >= 0, 0, 1)
  )

# x-axis limits with padding
x_min <- min(nnh_df$CI_low_pct, 0) - 1.5
x_max <- max(nnh_df$CI_high_pct, 0) + 3

# Add the "Assistant" font
font_add_google("Assistant", "Assistant")

# Automatically use showtext for new plots
showtext_auto()
# ----------------------------
# 3) Plot
# ----------------------------
gg <- ggplot(nnh_df, aes(x = RD_pct, y = pair)) +
  # Null (no-effect) zone: ±1% band
  geom_rect(
    aes(xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf),
    fill = "#F2F2F2",
    inherit.aes = FALSE
  ) +
  # Zero line
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.6) +

  # CIs
  geom_errorbarh(
    aes(xmin = CI_low_pct, xmax = CI_high_pct, alpha = alpha_val),
    height = 0.18,
    color = "black",
    linewidth = 0.9,
    show.legend = FALSE
  ) +

  # Points (colored by outcome group; dimmed if CI spans zero)
  geom_point(
    aes(color = group, alpha = alpha_val),
    size = 4.2,
    stroke = 0,
    show.legend = FALSE
  ) +

  # Labels (NNH/NNB or "No clear effect")
  geom_text(
    aes(x = label_x, label = NNH_label, hjust = label_hj),
    size = 8,
    color = "grey10"
  ) +

  scale_color_manual(
    values = c(
      "ESBL resistance" = pal_esbl,
      "Nitrofurantoin resistance" = pal_nitro
    )
  ) +
  scale_x_continuous(
    name = "Risk difference for resistance (%)  (Matched vs Reference)",
    limits = c(x_min, x_max),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  ylab(NULL) +
  facet_wrap(~group, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 24) +
  theme(
    text = element_text(family = "Assistant"), # Set the font family
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 12))
  )

print(gg)

# ----------------------------
# 4) Save (poster-ready)
# ----------------------------
ggsave(
  "admin/poster/ESBL_vs_Nitro_RD_NNH_poster.svg",
  gg,
  width = 12,
  height = 9,
  dpi = 320,
)
ggsave(
  "admin/poster/ESBL_vs_Nitro_RD_NNH_poster.png",
  gg,
  width = 12,
  height = 9,
  dpi = 320,
)

ggsave(
  "admin/poster/ESBL_vs_Nitro_RD_NNH_poster.jpeg",
  gg,
  width = 10,
  height = 6,
  dpi = 320,
)
