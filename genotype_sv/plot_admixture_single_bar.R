library(ggplot2)
library(reshape2)
library(dplyr)

# ================
# 参数：修改这里就行
# ================
K <- 6
Qfile <- sprintf("XX:/XX/sv_plink_pruned.%d.Q", K)
famfile <- "XX:/XX/admixture/sv_plink_pruned.fam"

# 读取
Q <- read.table(Qfile, stringsAsFactors = FALSE)
fam <- read.table(famfile, stringsAsFactors = FALSE)
sample_ids <- fam$V2
colnames(Q) <- paste0("Anc", 1:ncol(Q))
Q$SampleID <- sample_ids
Q$Group <- sub("_.*", "", Q$SampleID)

# 颜色分配（不够可自己扩展）
ancestry_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#4C9BEF",
  "#812f33", "#66a61e", "#e6ab02", "#9FD1EE"
)[1:K]

# 按 group 排序，group 内按 sample 排序
Q <- Q %>% arrange(Group, SampleID)

# 重塑
Qm <- melt(Q, id.vars=c("SampleID", "Group"),
           measure.vars=paste0("Anc",1:K))

# 个体顺序
Qm$SampleOrder <- factor(Qm$SampleID, levels=unique(Q$SampleID))

# ====================
# 生成条形图（单条带）
# ====================
p <- ggplot(Qm, aes(x=SampleOrder, y=value, fill=variable)) +
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = ancestry_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top"
  ) +
  labs(x="", y="Ancestry proportion", fill="Cluster")

# =========================
# 添加群体分隔线 & 群体标签
# =========================

# 计算每个 group 的区间
group_info <- Q %>% group_by(Group) %>% summarise(n = n())
group_info$cum_n <- cumsum(group_info$n)
group_info$start <- group_info$cum_n - group_info$n + 1
group_info$mid <- (group_info$start + group_info$cum_n)/2

# 加 group 分隔竖线
for (x in group_info$cum_n) {
  p <- p + geom_vline(xintercept=x+0.5, color="black", size=0.3)
}

# 添加 group 名称（侧边文本）
p <- p + annotate("text",
                  x = group_info$mid,
                  y = 1.05,
                  label = group_info$Group,
                  angle = 90,
                  size = 3)

# =======================
# 输出 PDF 和 PNG
# =======================
pdf(sprintf("admixture_K%d_single_bar.pdf", K), width=14, height=1.8)
print(p)
dev.off()

ggsave(sprintf("admixture_K%d_single_bar.png", K),
       plot=p, width=14, height=1.8, dpi=300)
