library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)
library(gtools)

# 文件路径
main_file <- "XX:/XX/DI_SV_siteLevel_chr.sorted.tsv"
highlight_file <- "XX:/XX/Manhattan_highlight_sites.tsv"

# 读入主数据
df <- fread(main_file, header = TRUE, sep = "\t")
colnames(df)[1:5] <- c("CHR", "POS", "AF1", "AF2", "DI")

# 去掉"chr"前缀（保持数字型排序）
df$CHR <- gsub("chr", "", df$CHR)
df$CHR <- factor(df$CHR, levels = mixedsort(unique(df$CHR)))

# 过滤掉DI缺失
df <- df[!is.na(df$DI), ]

# =======================
# Step 2. 累积基因组坐标 (修复整数溢出)
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
threshold <- quantile(df$DI, 0.95, na.rm = TRUE)

# =======================
# Step 4. 读取标注位点 (修复数据类型问题)
# =======================
# 正确读取highlight文件（2列，用空格分隔）
highlight_sites <- fread(highlight_file, header = FALSE, sep = " ")
cat("highlight文件列数:", ncol(highlight_sites), "\n")
print(head(highlight_sites))

# 如果列名不正确，手动设置
if (ncol(highlight_sites) == 2) {
  highlight_sites$V1 <- gsub("chr", "", highlight_sites$V1)
  highlight_sites$V2 <- as.numeric(highlight_sites$V2)  # 确保POS是数值型
  colnames(highlight_sites) <- c("CHR", "POS")
} else {
  # 如果列数不对，尝试其他分隔符
  highlight_sites <- fread(highlight_file, header = FALSE)
  cat("重新读取后的列数:", ncol(highlight_sites), "\n")
  if (ncol(highlight_sites) >= 2) {
    highlight_sites <- highlight_sites[, 1:2]
    highlight_sites$V1 <- gsub("chr", "", highlight_sites$V1)
    highlight_sites$V2 <- as.numeric(highlight_sites$V2)
    colnames(highlight_sites) <- c("CHR", "POS")
  } else {
    stop("无法正确读取highlight文件")
  }
}

# 确保数据类型匹配
df$POS <- as.numeric(df$POS)
highlight_sites$POS <- as.numeric(highlight_sites$POS)

# 取出匹配位点
df_highlight <- df %>%
  inner_join(highlight_sites, by = c("CHR", "POS"))

cat("读取到标注位点：", nrow(df_highlight), "个\n")

# =======================
# Step 5. 绘制曼哈顿图
# =======================
p <- ggplot(df, aes(x = BP_cum, y = DI, color = as.factor(as.numeric(CHR) %% 2))) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("#4D96C2", "#9DBCDA")) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "#3498db", size = 0.8) +
  geom_point(data = df_highlight, aes(x = BP_cum, y = DI),
             shape = 24, color = "#bf0b00", fill = NA, size = 3, stroke = 1.2) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  labs(x = "Chromosome", y = "DI value",
       title = "Manhattan Plot of SV-based DI (Top 5% dashed)") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# =======================
# Step 6. 输出PDF
# =======================
pdf("DI_SV_Manhattan_highlight_yd.pdf", width = 11, height = 5)
print(p)
dev.off()

cat("DI_SV_Manhattan_highlight.pdf\n")