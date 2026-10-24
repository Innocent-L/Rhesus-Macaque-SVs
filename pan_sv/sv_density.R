library(data.table)
library(ggplot2)

sizes <- fread("XX/XX/genome.sizes",
               col.names = c("chr", "length"))
svdens <- fread("XX/XX/SV_density_1Mb.bed",
                col.names = c("chr", "start", "end", "count"))

chr_order <- paste0("chr", 1:21)
chr_order_rev <- rev(chr_order)

sizes <- sizes[chr %in% chr_order]
svdens <- svdens[chr %in% chr_order]

sizes[, chr := factor(chr, levels = chr_order_rev)]
svdens[, chr := factor(chr, levels = chr_order_rev)]

chr_bars <- sizes[, .(chr, start = 0, end = length)]

p <- ggplot() +
  geom_rect(
    data = svdens,
    aes(
      xmin = as.numeric(chr) - 0.25,
      xmax = as.numeric(chr) + 0.25,
      ymin = start,
      ymax = end,
      fill = count
    ),
    color = NA
  ) +
  geom_rect(
    data = chr_bars,
    aes(
      xmin = as.numeric(chr) - 0.25,
      xmax = as.numeric(chr) + 0.25,
      ymin = start,
      ymax = end
    ),
    fill = NA,
    color = "grey30",
    linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = c("#fff7ec", "#fdd49e", "#fdbb84", "#ef6548", "#b30000","#470000"),
    name = "SV count"
  ) +
  scale_x_continuous(
    breaks = 1:length(chr_order_rev),
    labels = chr_order_rev,
    expand = c(0.02, 0.02)
  ) +
  coord_flip() +
  labs(title = "Genome-wide SV density in 1-Mb windows") +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("XX/XX/SV_density_1Mb_publication.pdf",
       p, width = 10, height = 6)