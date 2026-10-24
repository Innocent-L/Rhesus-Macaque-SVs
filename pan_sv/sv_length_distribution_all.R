library(dplyr)
library(ggplot2)

# 读取数据
sv <- read.table("X:/XX/sv_lengths.tsv", header = TRUE, sep = "\t")

# 主图：0-2kb分布
# 转为 kb
sv$length_kb <- sv$LENGTH / 1000

# 设置分组区间：0-2kb 每 0.1kb，一个 bin；>2kb 为最后一组
breaks <- c(seq(0, 2, by = 0.1), Inf)
labels <- c(
  paste0(seq(0, 1.9, 0.1), "-", seq(0.1, 2, 0.1), "kb"),
  ">2kb"
)

sv$bin <- cut(sv$length_kb, breaks = breaks, labels = labels, right = FALSE)

# 统计每个 bin 数量
sv_bin_count <- sv %>% group_by(bin) %>% summarise(count = n())

# 输出主图统计表格
write.table(sv_bin_count, "sv_bin_count.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("主图统计表格已保存为: sv_bin_count.tsv\n")

# 绘制主图
p_main <- ggplot(sv_bin_count, aes(x = bin, y = count)) +
  geom_bar(stat = "identity", fill = "grey30") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 5)
  ) +
  labs(
    x = "SV size (kb)",
    y = "Count",
    title = "SV Length Distribution (0–2kb, >2kb merged)"
  )

# 保存主图PDF
ggsave("SV_length_distribution_main.pdf", p_main, width = 3, height = 5)
cat("主图已保存为: SV_length_distribution_main.pdf\n")

# 插图：300-400 bp 精细分布
# 300-400 bp 区间提取
sv_300_400 <- sv %>% 
  filter(LENGTH >= 300 & LENGTH <= 400)

# 输出300-400bp的原始数据
write.table(sv_300_400, "sv_300_400_raw.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("300-400bp原始数据已保存为: sv_300_400_raw.tsv\n")

# 按 1 bp 统计数量
sv_300_400_count <- sv_300_400 %>%
  group_by(LENGTH) %>%
  summarise(count = n())

# 输出精细分布统计表格
write.table(sv_300_400_count, "sv_300_400_count.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("精细分布统计表格已保存为: sv_300_400_count.tsv\n")

# 绘制插图
p_inset <- ggplot(sv_300_400_count, aes(x = LENGTH, y = count)) +
  geom_line(color = "orange", size = 0.4) +
  theme_bw() +
  labs(
    x = "SV length (bp, 300–400)",
    y = "Count",
    title = "Fine-scale SV Length Distribution (300–400bp)"
  )

# 保存插图PDF
ggsave("SV_length_300_400_inset.pdf", p_inset, width = 2, height = 2)
cat("精细分布图已保存为: SV_length_300_400_inset.pdf\n")

# 输出数据摘要
cat("\n=== 数据摘要 ===\n")
cat("总SV数量:", nrow(sv), "\n")
cat("300-400bp区间SV数量:", nrow(sv_300_400), "\n")
cat("主图分组数量:", nrow(sv_bin_count), "\n")
cat("精细分布数据点数量:", nrow(sv_300_400_count), "\n")

# 保存完整的处理后的数据（可选）
write.table(sv, "sv_processed.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("完整处理后的数据已保存为: sv_processed.tsv\n")
