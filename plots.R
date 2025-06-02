# plot results


plt <- ggplot(res_df, aes(x = est, y = Estimator)) +
  geom_point(position = position_dodge(width = 0.5), size = 3.5) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                position = position_dodge(width = 0.5), width = 0.3, linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.5) +
  geom_vline(xintercept = 10, linetype = "dashed", linewidth = 1.5, color = "red") +
  labs(title = "ATE Estimates",
       x = "", y = "", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14))
plt