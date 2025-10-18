# Packages
# install.packages(c("ggplot2","ggtext","showtext","sysfonts"))
library(ggplot2)
library(ggtext)
library(showtext)
library(sysfonts)

# 1) Load the Dosis font (Google) once per session
use_dosis <- function() {
  try(font_add_google("Dosis", "Dosis"), silent = TRUE)
  showtext_auto() # render text with the added font
}
use_dosis()

# 2) Forest-plot function -----------------------------------------------
make_forest_plot <- function(
  plot_data,
  title = "",
  rd_x = NULL,
  x_limits = NULL,
  x_breaks = NULL,
  x_text = 22,
  family = "Dosis"
) {
  # Build colored HTML y-labels from the incoming data
  plot_data$y_label <- sprintf(
    "<span style='color:%s; font-weight:bold;'>%s</span> &#9656; %s",
    plot_data$hex_color,
    as.character(plot_data$group),
    plot_data$effect_on
  )

  # Set plot range
  if (is.null(x_limits)) {
    x_limits <- c(
      -5,
      max(plot_data$ci_high) + 5
    )
  }

  # Set x-axis breaks if not provided
  if (is.null(x_breaks)) {
    x_breaks <- seq(
      from = -5,
      to = ceiling(max(plot_data$ci_high) + 1),
      by = 5
    )
  }

  # Set RD text placement if not provided
  if (is.null(rd_x)) {
    rd_x <- max(plot_data$ci_high) + 1
  }

  # Keep the incoming order (top = first row)
  plot_data$y_label <- factor(
    plot_data$y_label,
    levels = rev(unique(plot_data$y_label))
  )

  p <- ggplot(plot_data, aes(x = point_estimate, y = y_label)) +
    # zero line
    geom_vline(xintercept = 0, color = "#343a40", linewidth = 0.5) +
    # CIs
    geom_linerange(
      aes(xmin = ci_low, xmax = ci_high, color = I(hex_color)),
      linewidth = 2,
      lineend = "round"
    ) +
    # Points
    geom_point(aes(color = I(hex_color)), size = 5) +
    # Right-side RD/NNH text
    geom_text(
      aes(label = rd_text),
      x = rd_x,
      hjust = 0,
      size = x_text - 14,
      lineheight = 0.9,
      family = family,
      color = "#343a40"
    ) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks, expand = c(0, 0)) +
    labs(x = "Risk Difference (%)", y = NULL) +
    theme_minimal(base_family = family) +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      plot.title = element_text(
        size = 30,
        face = "bold",
        hjust = 0.5,
        color = "#343a40",
        margin = margin(b = 25)
      ),
      axis.title.x = element_text(
        size = x_text + 2,
        color = "#343a40",
        margin = margin(t = 20)
      ),
      axis.text.x = element_text(size = x_text, color = "#343a40"),
      axis.text.y = element_markdown(
        family = family,
        size = x_text + 2,
        color = "#343a40",
        hjust = 1,
        lineheight = 1.1,
        margin = margin(r = -30)
      ),
      axis.line.x = element_line(color = "#343a40", linewidth = 0.5),
      axis.ticks.x = element_line(color = "#343a40"),
      axis.ticks.length.x = unit(0.25, "cm"),
      axis.ticks.y = element_blank()
    )

  return(p)
}

# --- Example usage (with your existing plot_data) ---
plot_data <- data.frame(
  # The order is reversed here to match the plot's top-to-bottom display
  group = factor(
    c("FQ", "AMC", "2nd-gen CEP"),
    levels = c("2nd-gen CEP", "AMC", "FQ")
  ),
  effect_on = c("1st-CEP", "1st-CEP", "1st-CEP"),
  point_estimate = c(22.1, 8.2, 7.8),
  ci_low = c(20, 6.1, 5.7),
  ci_high = c(24.2, 10.3, 9.7),
  hex_color = c("#e76ea2", "#f5a250", "#61b5cf"),
  # Text for the right-hand side annotations
  rd_text = c("RD 22.1%\n(NNH=4)", "RD 8.2%\n(NNH=12)", "RD 7.8%\n(NNH=12)")
)

p_1 <- make_forest_plot(plot_data)
print(p_1)
ggsave(
  "admin/ab_resist_presentation/trial1_forest.svg",
  p_1,
  bg = "transparent"
) # 16:9-ish

plot_data_2 <- data.frame(
  # The order is reversed here to match the plot's top-to-bottom display
  group = c("2nd-CEP", "2nd-CEP"),
  effect_on = c("FQ", "AMC"),
  point_estimate = c(-21, -5.1),
  ci_low = c(-23, -6.7),
  ci_high = c(-19.1, -3.4),
  hex_color = c("#61b5cf", "#61b5cf"),
  # Text for the right-hand side annotations
  rd_text = c("RD 21%\n(NNT=5)", "RD 5.1%\n(NNT=20)")
)

p_2 <- make_forest_plot(
  plot_data_2,
  rd_x = 1,
  x_limits = c(-26, 5),
  x_breaks = seq(-25, 5, 5)
)
print(p_2)
ggsave(
  "admin/ab_resist_presentation/trial2_forest.svg",
  p_2,
  width = 10,
  height = 5.625,
  units = "in",
  bg = "transparent"
)

plot_data_3 <- data.frame(
  # The order is reversed here to match the plot's top-to-bottom display
  group = c(
    "AMC",
    "2nd-CEP",
    "FQ",
    "2nd-CEP",
    "2nd-CEP"
  ),
  effect_on = c("1st-CEP", "1st-CEP", "1st-CEP", "AMC", "FQ"),
  point_estimate = c(3.7, 8.5, 5.3, -4.2, -2.1),
  ci_low = c(2, 6.8, 3.6, -5.6, -3.9),
  ci_high = c(5.3, 10.2, 7.1, -2.9, -0.4),
  hex_color = c("#f5a250", "#61b5cf", "#e76ea2", "#61b5cf", "#61b5cf"),
  # Text for the right-hand side annotations
  rd_text = c(
    "RD 3.7%\n(NNH=28)",
    "RD 8.5%\n(NNH=12)",
    "RD 5.3%\n(NNH=18)",
    "RD 4.2%\n(NNT=23)",
    "RD 2.1%\n(NNT=46)"
  )
)

p_3 <- make_forest_plot(
  plot_data_3,
  x_limits = c(-10, 15),
  x_breaks = seq(-10, 15, 5)
)
print(p_3)
ggsave(
  "admin/ab_resist_presentation/trial3_forest.svg",
  p_3,
  width = 10,
  height = 5.625,
  units = "in",
  bg = "transparent"
)
