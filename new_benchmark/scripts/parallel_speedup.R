library(ggplot2)

df <- read.table("results/parallel_speedup.tsv", header = TRUE, sep = "\t")

base_time <- df$elapsed_seconds[df$threads == min(df$threads)]
df$speedup <- base_time / df$elapsed_seconds

p <- ggplot(df, aes(x = threads, y = speedup)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Query parallel speedup",
    x = "Threads",
    y = "Speedup",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey10", linewidth = 0.2),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) +
  expand_limits(x = 0, y = 0)

ggsave("results/parallel_speedup.pdf", p, width = 5, height = 5, dpi = 100, bg = "white")
print(p)
