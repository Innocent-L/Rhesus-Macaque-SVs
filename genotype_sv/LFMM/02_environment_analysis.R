# 02: independent-site geographic/BIO correlations and climate PCA
# Set this path to the location of site_bio.xlsx on your Windows computer.
SITE_BIO <- "XX:/XX/site_bio.xlsx"
OUT <- "results_02_environment"
packages <- c("readxl", "ggplot2")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Install packages: install.packages(c(", paste(sprintf('"%s"',missing), collapse=","), "))")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
s <- as.data.frame(readxl::read_excel(SITE_BIO)); names(s) <- tolower(trimws(names(s)))
bios <- paste0("bio",1:19); req <- c("site_id","lat","lon",bios)
if (!all(req %in% names(s))) stop("Missing site_bio columns: ", paste(setdiff(req,names(s)),collapse=","))
if (anyNA(s[,req])) stop("Missing site/coordinate/BIO values; inspect site_bio.xlsx")
if (anyDuplicated(s$site_id)) stop("Duplicate site_id rows; do not silently discard")
for (nm in c("lat","lon",bios)) s[[nm]] <- as.numeric(s[[nm]])
if (anyNA(s[,c("lat","lon",bios)])) stop("Non-numeric coordinate/BIO value")
rownames(s) <- s$site_id
cat("Number of independent sites:",nrow(s),"\n")
write.table(s[,req], file.path(OUT,"02_sites_used.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
X <- s[,c("lat","lon",bios)]
for (method in c("pearson","spearman")) {
  C <- cor(X, method=method)
  write.csv(C, file.path(OUT,paste0("02_geo_BIO_",method,"_cor.csv")))
  pdf(file.path(OUT,paste0("02_geo_BIO_",method,"_heatmap.pdf")),width=11,height=10)
  heatmap(C, symm=TRUE, margins=c(7,7), main=paste("11-site geography + BIO",method,"correlations"))
  dev.off()
}
lat_bio <- do.call(rbind, lapply(bios,function(nm){
  a <- suppressWarnings(cor.test(s$lat,s[[nm]],method="pearson"))
  b <- suppressWarnings(cor.test(s$lat,s[[nm]],method="spearman",exact=FALSE))
  data.frame(BIO=nm,pearson_r=unname(a$estimate),pearson_p=a$p.value,
             spearman_rho=unname(b$estimate),spearman_p=b$p.value)
}))
lat_bio$pearson_BH <- p.adjust(lat_bio$pearson_p,"BH")
lat_bio$spearman_BH <- p.adjust(lat_bio$spearman_p,"BH")
lat_bio <- lat_bio[order(-abs(lat_bio$pearson_r)),]
write.table(lat_bio,file.path(OUT,"02_latitude_BIO_correlations.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
if (any(apply(s[,bios],2,sd)==0)) stop("Some BIO variables are constant across sites")
pca <- prcomp(s[,bios],center=TRUE,scale.=TRUE)
variance <- pca$sdev^2/sum(pca$sdev^2)
write.table(data.frame(PC=paste0("PC",seq_along(variance)),variance=variance,cumulative=cumsum(variance)),
 file.path(OUT,"02_climate_PCA_variance.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(data.frame(BIO=rownames(pca$rotation),pca$rotation,check.names=FALSE),
 file.path(OUT,"02_climate_PCA_loadings.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
scores <- data.frame(site_id=s$site_id,pca$x,lat=s$lat,lon=s$lon,check.names=FALSE)
write.table(scores,file.path(OUT,"02_climate_PCA_site_scores.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
pdf(file.path(OUT,"02_climate_PCA_sites.pdf"),width=7,height=6)
plot(scores$PC1,scores$PC2,pch=19,xlab=sprintf("Climate PC1 (%.1f%%)",variance[1]*100),
 ylab=sprintf("Climate PC2 (%.1f%%)",variance[2]*100))
text(scores$PC1,scores$PC2,labels=scores$site_id,pos=3,cex=.8)
dev.off()
write.table(data.frame(PC=c("PC1","PC2"),r_lat=c(cor(scores$PC1,s$lat),cor(scores$PC2,s$lat)),
 r_lon=c(cor(scores$PC1,s$lon),cor(scores$PC2,s$lon))),
 file.path(OUT,"02_climatePC_geography_cor.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
for (nm in head(lat_bio$BIO,2)) {
 p <- ggplot2::ggplot(s,ggplot2::aes(x=.data$lat,y=.data[[nm]])) +
 ggplot2::geom_point(size=2.7)+ggplot2::geom_text(ggplot2::aes(label=.data$site_id),vjust=-0.6,check_overlap=TRUE)+
 ggplot2::geom_smooth(method="lm",se=TRUE)+ggplot2::theme_classic()+
 ggplot2::labs(x="Latitude (degrees)",y=nm,title=paste("Latitude vs",nm))
 ggplot2::ggsave(file.path(OUT,paste0("02_lat_vs_",nm,".pdf")),p,width=6,height=5)
}
cat("Environment outputs in", OUT, "\n")

#去除一个个体的位点为后续RDA分析提供PCA结果
library(readxl)

site <- read_excel("XX:/XX/site_bio.xlsx")

site10 <- subset(site, site_id != "zm4")

bio_cols <- paste0("bio", 1:19)

pca10 <- prcomp(
  site10[, bio_cols],
  center = TRUE,
  scale. = TRUE
)

scores10 <- as.data.frame(pca10$x)
scores10$site_id <- site10$site_id
scores10$lat <- site10$lat
scores10$lon <- site10$lon

write.table(
  scores10,
  "02_climate_PCA_site_scores_10sites.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    BIO = rownames(pca10$rotation),
    pca10$rotation
  ),
  "02_climate_PCA_loadings_10sites.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
