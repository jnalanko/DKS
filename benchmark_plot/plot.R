# Load required library
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)

# Read the TSV file (update the filename as needed)
data <- read.table("results.csv", header = TRUE, sep = ",")

# Ensure proper data types
data$k <- as.factor(data$k)
data$threads <- as.integer(data$threads)
data$time <- as.numeric(data$time) / 60
data$mem <- as.numeric(data$mem)
data$disk <- as.numeric(data$disk)
data$mode <- as.factor(data$mode)
data$query_time_per_base_ns <- as.numeric(data$query_time_per_base_ns)
data$mode_and_k <- interaction(data$mode, data$k, sep = "-k")
data$query_throughput <- 1 / data$query_time_per_base_ns * 1e9 # Bases / s

# --- Plot 1: Parallel Speedup vs Threads ---
p <- ggplot(data, aes(x = threads, y = time, shape = mode, color = k)) +
  geom_point(size = 2) +
  geom_line() +
  guides(color = guide_legend(override.aes = list(shape = NA))) +
  labs(
    title = "Indexing time",
    x = "Threads",
    y = "Time (minutes)",
  ) +
  theme_minimal(base_size = 14) +
  scale_shape_manual(
    name = "Use disk",                     # legend title
    values = c("disk" = 16, "mem" = 17),                 # shape codes
    labels = c("Yes", "No")      # legend labels
  ) +
  theme(
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey10", linewidth = 0.2),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5), # Horitonzal centering
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) +
  expand_limits(x = 0, y = 0)

# Save the plot
ggsave("time_vs_threads.png", p, width = 5, height = 5, dpi = 100, bg = "white")
print(p)

custom_titles <- c(
  "disk" = "Disk-based construction",
  "mem" = "In-memory construction"
)

# --- Compute speedup: baseline = time at MIN(threads) per (k, mode) ---
baseline <- data %>%
  group_by(k, mode) %>%
  filter(threads == min(threads)) %>%
  summarise(base_time = median(time, na.rm = TRUE), .groups = "drop")

data_speedup <- data %>%
  left_join(baseline, by = c("k", "mode")) %>%
  mutate(speedup = base_time / time)

# --- Plot 3: Parallel Speedup vs Threads ---
# Ideal linear speedup is y = x (dashed reference line)
p_speedup <- ggplot(data_speedup, aes(x = threads, y = speedup, shape = mode, color = k)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line() + 
  geom_point(size = 2) +
  guides(color = guide_legend(override.aes = list(shape = NA))) +
  labs(
    title = "Indexing parallel speedup",
    x = "Threads",
    y = "Speedup",
  ) +
  scale_shape_manual(
    name = "Use disk",                     # legend title
    values = c("disk" = 16, "mem" = 17),                 # shape codes
    labels = c("Yes", "No")      # legend labels
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey10", linewidth = 0.2),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) +
  expand_limits(x = 32, y = 32)

ggsave("benchmark_speedup.png", p_speedup, width = 5, height = 5, dpi = 100, bg = "white")
print(p_speedup)


# --- Plot 4: Time-space plane
timespace_threads = 8
p_timespace <- ggplot(data %>% filter(threads == timespace_threads), aes(x = time, y = mem / 2^30, color = k, shape = mode)) +
  geom_point(size = 4) +
  labs(
    title = paste0("Indexing time vs memory (",  timespace_threads, " threads)"),
    x = "Time (min)",
    y = "Memory (GiB)",
    shape = "Mode",
    color = "k"
  ) +
  scale_shape_manual(
    name = "Use disk",                     # legend title
    values = c("disk" = 16, "mem" = 17),                 # shape codes
    labels = c("Yes", "No")      # legend labels
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

print(p_timespace)
ggsave("time_mem.png", p_timespace, width = 5, height = 5, dpi = 100, bg = "white")

# -- Plot 5: Query Throughput
p_query <- ggplot(data %>% filter(mode == "mem"), aes(x = threads, y = query_throughput / 1e6, color = k)) +
  geom_line() +
  geom_abline(slope = min((data %>% filter(k == 31))$query_throughput) / 1e6, intercept = 0, linetype = "dashed", color = "#BF4040") + # Ideal speedup
  geom_abline(slope = min((data %>% filter(k == 63))$query_throughput) / 1e6, intercept = 0, linetype = "dashed", color = "#00BFC4") + # Ideal speedup
  geom_point(size = 2) +
  labs(
    title = paste0("Query throughput"),
    x = "Threads",
    y = "Throughput (Mbases / s)",
    color = "k"
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

print(p_query)

p_grid = grid.arrange(p, p_speedup, p_timespace, p_query, ncol = 2, top = textGrob("Human genome", gp = gpar(fontsize = 24, fontface = "bold")))
print(p_grid)
ggsave("benchmarks_combined.png", p_grid, width = 9, height = 9, dpi = 200, bg = "white")