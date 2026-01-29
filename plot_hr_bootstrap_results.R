library(ggplot2)
library(dplyr)

# Load HR results from hr_point_est_with_bootstrap_ci.R
hr_file <- "hr_point_est_with_bootstrap_ci.csv"

if (!file.exists(hr_file)) {
  stop(sprintf("HR results file not found: %s\nPlease run hr_point_est_with_bootstrap_ci.R first.", hr_file))
}

hr_df <- read.csv(hr_file)

# Map model names to cleaner names matching Censoring_RWD_HR.R style
# Define the desired order (top to bottom after coord_flip)
desired_order <- c(
  "CUMC LSPS-IPTW and IPCW",
  "CUMC IPCW",
  "CUMC LSPS-IPTW",
  "CUMC Unadjusted"
)

hr_plot_df <- hr_df %>%
  mutate(
    names = case_when(
      Model == "unadjusted" ~ "CUMC Unadjusted",
      Model == "iptw_trunc" ~ "CUMC LSPS-IPTW",
      Model == "ipcw_stab_trunc" ~ "CUMC IPCW",
      Model == "combined" ~ "CUMC LSPS-IPTW and IPCW",
      TRUE ~ Model
    ),
    hr = HR_estimate,
    lower = HR_CI_lower,
    upper = HR_CI_upper
  ) %>%
  select(names, hr, lower, upper)

# Optional: Add legend row (uncomment and adjust if needed)
# Legend will appear at the very top if uncommented
legend_row <- data.frame(
  names = "LEGEND LSPS-Matching",
  hr = 0.84,  # Update with actual legend value if needed
  lower = 0.75,  # Update with actual legend CI if needed
  upper = 0.95   # Update with actual legend CI if needed
)
hr_plot_df <- rbind(legend_row, hr_plot_df)
desired_order <- c(desired_order, "LEGEND LSPS-Matching")

# Set factor levels to desired order (top to bottom after coord_flip)
hr_plot_df$names <- factor(hr_plot_df$names, levels = desired_order)

# Define colors for each bar
# If legend row is included, make it red; otherwise all black
colors <- ifelse(hr_plot_df$names == "LEGEND LSPS-Matching", "red", "black")

# Create the forest plot
fp <- ggplot(data = hr_plot_df, aes(x = names, y = hr, ymin = lower, ymax = upper)) +
  geom_pointrange(aes(color = colors)) + 
  geom_hline(yintercept = 1, lty = 2) +  # Add a dotted line at HR=1
  coord_flip() +  # Flip coordinates (puts labels on y axis)
  xlab("Method") + 
  ylab("Hazard Ratio") +
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
ggsave("results/HRplots.png", fp, 
       width = 14,
       height = 7,
       dpi = 600)

ggsave("results/HRplots.pdf", fp, 
       width = 14, 
       height = 7,
       dpi = 600)

# 
# # FOR PAPER
# ggsave("results/HRplots.png", fp, 
#        width = 14,
#        height = 7,
#        dpi = 600)
# 
# ggsave("results/HRplots.png", fp, 
#        width = 10,
#        height = 6,
#        dpi = 600)

cat("\nHR forest plot saved successfully!\n")
