library(ggplot2)
library(dplyr)

# Load RMST results from rmst_point_est_with_bootstrap_ci.R
rmst_file <- "rmst_point_est_with_bootstrap_ci.csv"

if (!file.exists(rmst_file)) {
  stop(sprintf("RMST results file not found: %s\nPlease run rmst_point_est_with_bootstrap_ci.R first.", rmst_file))
}

rmst_df <- read.csv(rmst_file)

# Get tau value from the dataframe (should be the same for all rows)
tau <- unique(rmst_df$tau)[1]

# Map model names to cleaner names matching Censoring_RWD_HR.R style
# Define the desired order (top to bottom after coord_flip)
desired_order <- c(
  "CUMC LSPS-IPTW and IPCW",
  "CUMC IPCW",
  "CUMC LSPS-IPTW",
  "CUMC Unadjusted"
)

rmst_plot_df <- rmst_df %>%
  mutate(
    names = case_when(
      Model == "unadjusted" ~ "CUMC Unadjusted",
      Model == "iptw" ~ "CUMC LSPS-IPTW",
      Model == "ipcw" ~ "CUMC IPCW",
      Model == "combined" ~ "CUMC LSPS-IPTW and IPCW",
      TRUE ~ Model
    ),
    estimate = RMST_diff_estimate,
    lower = RMST_diff_CI_lower,
    upper = RMST_diff_CI_upper
  ) %>%
  select(names, estimate, lower, upper)

# Optional: Add legend row (uncomment and adjust if needed)
# legend_row <- data.frame(
#   names = "LEGEND LSPS-Matching",
#   estimate = 0.0,  # Update with actual legend value if needed
#   lower = -10.0,  # Update with actual legend CI if needed
#   upper = 10.0   # Update with actual legend CI if needed
# )
# rmst_plot_df <- rbind(legend_row, rmst_plot_df)
# desired_order <- c("LEGEND LSPS-Matching", desired_order)

# Set factor levels to desired order (top to bottom after coord_flip)
rmst_plot_df$names <- factor(rmst_plot_df$names, levels = desired_order)

# Define colors for each bar
# If legend row is included, make it red; otherwise all black
colors <- ifelse(rmst_plot_df$names == "LEGEND LSPS-Matching", "red", "black")

# Create the forest plot
fp <- ggplot(data = rmst_plot_df, aes(x = names, y = estimate, ymin = lower, ymax = upper)) +
  geom_pointrange(aes(color = colors)) + 
  geom_hline(yintercept = 0, lty = 2) +  # Add a dotted line at RMST difference=0
  coord_flip() +  # Flip coordinates (puts labels on y axis)
  xlab("Method") + 
  ylab(sprintf("RMST Difference (tau = %d)", tau)) +
  scale_color_identity() +  # Use the specified colors directly
  theme_bw() +  # Use a white background
  theme(
    plot.title = element_text(size = 20),
    axis.text = element_text(size = 26),
    axis.title = element_text(size = 26)
  )

print(fp)

# Save plots
# FOR AMIA SUBMISSION
ggsave("results/RMSTplots.png", fp, 
       width = 9,
       height = 5,
       dpi = 600)

ggsave("results/RMSTplots.pdf", fp, 
       width = 14,
       height = 7,
       dpi = 600)


# # FOR PAPER
# ggsave("results/RMSTplots.png", fp, 
#        width = 14,
#        height = 7,
#        dpi = 600)
# 
# ggsave("results/RWD_RMSTplots.pdf", fp, 
#        width = 10,
#        height = 6,
#        dpi = 600)

cat("\nRMST forest plot saved successfully!\n")
