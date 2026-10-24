# 加载必要包
library(ggplot2)
library(plotly)
library(dplyr)

# ----------------------------------------------
# 1. 读取 PCA eigenvec
# ----------------------------------------------
eigenvec_file <- "XX:/XX/sv_pca_maf005_pruned.eigenvec"
eigenval_file <- "XX:/XX/sv_pca_maf005_pruned.eigenval"

pca_data <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE)

# 自动命名列，第一二列是 sample 信息
pc_count <- ncol(pca_data) - 2
colnames(pca_data) <- c("SampleID", "SampleID2", paste0("PC", 1:pc_count))

# 群体名提取（样本ID '_' 前的部分）
pca_data$Group <- sub("_.*", "", pca_data$SampleID)


# ----------------------------------------------
# 2. 读取 eigenval（用于计算方差贡献率）
# ----------------------------------------------
eigenval <- scan(eigenval_file)

# 每个 PC 的解释度
var_explained <- eigenval / sum(eigenval)

pc1_percent <- round(var_explained[1] * 100, 2)
pc2_percent <- round(var_explained[2] * 100, 2)
pc3_percent <- round(var_explained[3] * 100, 2)

cat("PC1 explained:", pc1_percent, "%\n")
cat("PC2 explained:", pc2_percent, "%\n")
cat("PC3 explained:", pc3_percent, "%\n")


# ----------------------------------------------
# 3. 自定义颜色
# ----------------------------------------------
group_colors <- c(
  "cx" = "#f8766d",
  "fj" = "#c5b233",
  "hb" = "#33c860",
  "hn" = "#33ccd0",
  "yd" = "#629dff",
  "zm" = "#bc51ff"
)


# ----------------------------------------------
# 4. 绘制 2D PCA 图（带百分比）
# ----------------------------------------------
p2d <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  theme_classic() +
  labs(
    title = "PCA 2D Plot",
    x = paste0("PC1 (", pc1_percent, "%)"),
    y = paste0("PC2 (", pc2_percent, "%)")
  ) +
  theme(text = element_text(size = 14))

ggsave("PCA_2D_with_variance_maf005_pruned.pdf", p2d, width = 6, height = 5)


# ----------------------------------------------
# 5. 绘制 3D PCA 图（带百分比）
# ----------------------------------------------
if (pc_count >= 3) {
  p3d <- plot_ly(
    pca_data,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color = ~Group,
    colors = group_colors,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 4)
  ) %>%
    layout(
      title = "PCA 3D Plot",
      scene = list(
        xaxis = list(title = paste0("PC1 (", pc1_percent, "%)")),
        yaxis = list(title = paste0("PC2 (", pc2_percent, "%)")),
        zaxis = list(title = paste0("PC3 (", pc3_percent, "%)"))
      )
    )
  
  htmlwidgets::saveWidget(p3d, "PCA_3D_with_variance_maf005_pruned.html")
} else {
  cat("⚠️ 你的 PCA 文件没有 PC3，无法绘制 3D 图\n")
}


# ----------------------------------------------
# 显示图
# ----------------------------------------------
print(p2d)
if (exists("p3d")) p3d
# ----------------------------------------------
# 网页调整角度后输出pdf图
#F12后打开控制台，输入document.querySelectorAll('div.js-plotly-plot')[0].layout.scene.camera.eye
#得到X\Y\Z参数.

# ----------------------------------------------


