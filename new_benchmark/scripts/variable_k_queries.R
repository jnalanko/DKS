library(ggplot2)

df <- read.table("results/variable_k_queries.tsv", header = TRUE, sep = "\t")

p <- ggplot(df, aes(x = k, y = query_time_ns_per_bp)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Single-threaded query time",
    x = "k",
    y = "ns / character",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey10", linewidth = 0.2),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) +
  expand_limits(y = 0)

ggsave("results/query_time_vs_k.pdf", p, width = 5, height = 5, dpi = 100, bg = "white")
print(p)
