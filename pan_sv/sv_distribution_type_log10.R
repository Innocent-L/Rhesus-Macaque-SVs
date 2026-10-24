library(readr)
library(dplyr)
library(ggplot2)

# 读取数据
sv <- read_tsv("XX:/XX/sv_lengths.tsv")

# 自定义颜色
my_colors <- c(
  "DEL"   = "#73ABCF",
  "INS"   = "#ED746A",
  "INV"   = "#FAAE5F",
  "DUP"   = "#354898",
  "CONTR" = "#972D36"
)

# 数据整理：直接取 log10
sv_plot <- sv %>%
  filter(
    !is.na(TYPE),
    !is.na(LENGTH),
    LENGTH > 0,
    TYPE %in% names(my_colors)
  ) %>%
  mutate(
    TYPE = factor(TYPE, levels = c("INS", "DEL", "DUP", "INV", "CONTR")),
    log10_LENGTH = log10(LENGTH)
  )

# 作图
p <- ggplot(sv_plot, aes(x = TYPE, y = log10_LENGTH, fill = TYPE)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.65,
    color = "black",
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.size = 0.25,
    outlier.alpha = 0.35,
    color = "black",
    linewidth = 0.35,
    alpha = 0.9
  ) +
  scale_fill_manual(values = my_colors) +
  labs(
    x = "SV type",
    y = expression(log[10]~"(SV length, bp)"),
    title = "SV length distribution"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(8, 8, 8, 8)
  )

p

ggsave(
  "SV_length_distribution_log10_real_values.pdf",
  p,
  width = 8,
  height = 6
)