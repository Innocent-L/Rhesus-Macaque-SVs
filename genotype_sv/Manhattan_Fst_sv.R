library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(gtools)
library(tidyr)

# 文件路径
main_file <- "XX:/XX/SV_Fst.weir_chr.fst"
highlight_file <- "XX:/XX/Manhattan_highlight_sites.tsv"

# 读入主数据 - 现在只有3列
df <- fread(main_file, header = FALSE, sep = "\t")
colnames(df)[1:3] <- c("CHR", "POS", "FST")

# 去掉"chr"前缀（保持数字型排序）
df$CHR <- gsub("chr", "", df$CHR)
df$CHR <- factor(df$CHR, levels = mixedsort(unique(df$CHR)))

# 将负FST值设为0
df$FST_nonneg <- ifelse(df$FST < 0, 0, df$FST)

# 过滤掉FST缺失值
df <- df[!is.na(df$FST_nonneg), ]

# 打印FST值的统计信息
cat("原始FST值统计:\n")
cat("最小值:", min(df$FST, na.rm = TRUE), "\n")
cat("最大值:", max(df$FST, na.rm = TRUE), "\n")
cat("均值:", mean(df$FST, na.rm = TRUE), "\n")
cat("负值数量:", sum(df$FST < 0, na.rm = TRUE), "\n\n")

cat("非负FST值统计:\n")
cat("最小值:", min(df$FST_nonneg, na.rm = TRUE), "\n")
cat("最大值:", max(df$FST_nonneg, na.rm = TRUE), "\n")
cat("均值:", mean(df$FST_nonneg, na.rm = TRUE), "\n")
cat("95%分位数:", quantile(df$FST_nonneg, 0.95, na.rm = TRUE), "\n\n")

# =======================
# Step 2. 累积基因组坐标
# =======================
df <- df %>%
  group_by(CHR) %>%
  arrange(POS) %>%
  mutate(BP_cum = as.numeric(POS)) %>%
  ungroup()

# 染色体长度
chr_info <- df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(as.numeric(BP_cum))) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len))

df <- df %>%
  left_join(chr_info, by = "CHR") %>%
  mutate(BP_cum = as.numeric(POS) + tot)

# x轴标签位置
axis_df <- df %>%
  group_by(CHR) %>%
  summarise(center = (min(BP_cum) + max(BP_cum)) / 2)

# =======================
# Step 3. Top 5% 阈值线
# =======================
# 用于绘图：负FST显示为0
df$FST_nonneg <- ifelse(df$FST < 0, 0, df$FST)

# 用于计算阈值：删除NA和负值，与Python脚本保持一致
fst_for_cutoff <- df$FST[!is.na(df$FST) & df$FST >= 0]

threshold <- quantile(
  fst_for_cutoff,
  probs = 0.95,
  na.rm = TRUE,
  type = 7
)

cat("参与阈值计算的非负FST数量:", length(fst_for_cutoff), "\n")
cat("FST Top 5%阈值:", threshold, "\n")
# =======================
# Step 4. 读取标注位点
# =======================
# 读取highlight文件
highlight_sites <- fread(highlight_file, header = FALSE, sep = "\t")
cat("highlight文件原始列数:", ncol(highlight_sites), "\n")

# 如果只有一列，则分割
if (ncol(highlight_sites) == 1) {
  cat("检测到单列数据，正在分割为两列...\n")
  
  # 使用tidyr的separate函数分割列
  highlight_sites <- highlight_sites %>%
    separate(V1, into = c("CHR", "POS"), sep = "\\s+")  # 使用一个或多个空白字符分割
  
  cat("分割后的列数:", ncol(highlight_sites), "\n")
} else if (ncol(highlight_sites) >= 2) {
  # 如果有多列，取前两列
  highlight_sites <- highlight_sites[, 1:2]
  colnames(highlight_sites) <- c("CHR", "POS")
}

# 处理CHR列：去掉"chr"前缀
highlight_sites$CHR <- gsub("chr", "", highlight_sites$CHR)
# 确保POS是数值型
highlight_sites$POS <- as.numeric(highlight_sites$POS)

# 移除任何NA值
highlight_sites <- highlight_sites[!is.na(highlight_sites$POS), ]

cat("处理后的highlight数据前6行:\n")
print(head(highlight_sites))

# 确保主数据中的POS是数值型
df$POS <- as.numeric(df$POS)

# 取出匹配位点
df_highlight <- df %>%
  inner_join(highlight_sites, by = c("CHR", "POS"))

cat("读取到标注位点：", nrow(df_highlight), "个\n")

if (nrow(df_highlight) == 0) {
  cat("警告：没有找到匹配的标注位点\n")
  cat("highlight文件中的染色体示例：", head(unique(highlight_sites$CHR)), "\n")
  cat("主数据中的染色体示例：", head(levels(df$CHR)), "\n")
}

# =======================
# Step 5. 绘制曼哈顿图
# =======================
p <- ggplot(df, aes(x = BP_cum, y = FST_nonneg, color = as.factor(as.numeric(CHR) %% 2))) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("#b6e3cf", "#6ea48c")) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "#3b6150", size = 0.8) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  scale_y_continuous(limits = c(0, max(df$FST_nonneg, na.rm = TRUE) * 1.05)) +  # 确保y轴从0开始
  labs(x = "Chromosome", y = "Fst value",
       title = "Manhattan Plot of SV-based Fst (Top 5% dashed)") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# 如果有匹配的标注位点，添加它们
if (nrow(df_highlight) > 0) {
  p <- p + geom_point(data = df_highlight, aes(x = BP_cum, y = FST_nonneg),
                      shape = 24, color = "#bf0b00", fill = NA, size = 3, stroke = 1.2)
}

# =======================
# Step 6. 输出PDF
# =======================
pdf("Fst_SV_Manhattan_highlight_new.pdf", width = 11, height = 5)
print(p)
dev.off()

cat("曼哈顿图已保存为 FManhattan_highlight.pdf\n")