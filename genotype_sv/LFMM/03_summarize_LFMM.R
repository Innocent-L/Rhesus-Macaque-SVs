# 04: Summarize existing genome-wide LFMM output, without re-running LFMM.
LFMM_ALL <- "LFMM2_all_19bio_results.tsv"
K_ROBUST <- "LFMM2_K_robust_candidates.tsv"
K_STRICT <- "LFMM2_K_strict_robust_candidates.tsv"
OUT <- "results_04_LFMM_summary"
dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
a <- read.delim(LFMM_ALL,check.names=FALSE)
req <- c("analysis_sv_id","bio","q_global_19bio","q_within_bio")
if(!all(req %in% names(a)))stop("Missing LFMM columns: ",paste(setdiff(req,names(a)),collapse=","))
if (anyDuplicated(a[,c("analysis_sv_id","bio")])) stop("Duplicate SV x BIO rows. Check the LFMM ID mapping before summarizing")
read_k <- function(path){
 if(!file.exists(path))return(NULL)
 x<-read.delim(path,check.names=FALSE)
 if(!all(c("bio","analysis_sv_id") %in% names(x)))stop("K file missing bio/analysis_sv_id: ",path)
 unique(x[,c("bio","analysis_sv_id")])
}
k1<-read_k(K_ROBUST);k2<-read_k(K_STRICT)
count_sig<-function(df,test){
 length(unique(df$analysis_sv_id[which(test)]))
}
ss<-do.call(rbind,lapply(paste0("bio",1:19),function(b){
 x<-a[a$bio==b,,drop=FALSE]
 data.frame(bio=b,n_tested=count_sig(x,rep(TRUE,nrow(x))),
 tier1_global=count_sig(x,is.finite(x$q_global_19bio)&x$q_global_19bio<.05),
 within_bio_FDR=count_sig(x,is.finite(x$q_within_bio)&x$q_within_bio<.05),
 tier2_only=count_sig(x,is.finite(x$q_within_bio)&x$q_within_bio<.05 & x$q_global_19bio>=.05),
 exploratory_q05_10=count_sig(x,is.finite(x$q_within_bio)&x$q_within_bio>=.05 & x$q_within_bio<.10),
 K_robust=if(is.null(k1))NA_integer_ else count_sig(k1[k1$bio==b,,drop=FALSE],rep(TRUE,sum(k1$bio==b))),
 K_strict=if(is.null(k2))NA_integer_ else count_sig(k2[k2$bio==b,,drop=FALSE],rep(TRUE,sum(k2$bio==b))))
}))
write.table(ss,file.path(OUT,"04_BIO_association_counts.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
formal<-a[(is.finite(a$q_within_bio)&a$q_within_bio<.05)| (is.finite(a$q_global_19bio)&a$q_global_19bio<.05),]
write.table(data.frame(metric=c("unique_SVs_tested","unique_formal_SVs","formal_SV_BIO_pairs"),
 n=c(length(unique(a$analysis_sv_id)),length(unique(formal$analysis_sv_id)),nrow(formal))),
 file.path(OUT,"04_LFMM_global_summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
pdf(file.path(OUT,"04_BIO_formal_association_counts.pdf"),width=10,height=5)
barplot(ss$within_bio_FDR,names.arg=ss$bio,las=2,ylab="SVs at within-BIO FDR < 0.05",
 main="Detected SV associations by BIO (not variance explained)")
dev.off()
cat("LFMM summary complete. Do not interpret high counts as causal strength.\n")
