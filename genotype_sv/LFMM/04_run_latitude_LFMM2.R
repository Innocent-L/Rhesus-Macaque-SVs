# Latitude LFMM2, K=5 primary, K=4/6 sensitivity
# Install LEA if missing: BiocManager::install("LEA")
# Set INPUT_DIR to folder produced by 05a_prepare_latitude_LFMM.py
INPUT_DIR <- "XX:/XX/latitude_LFMM_input"
OUT_DIR <- "XX:/XX/results_05_latitude_LFMM"
K_VALUES <- c(4L,5L,6L)
PRIMARY_K <- 5L
set.seed(20261010)
if (!requireNamespace('LEA', quietly=TRUE)) stop('Please install LEA: BiocManager::install("LEA")')
dir.create(OUT_DIR,showWarnings=FALSE,recursive=TRUE)
read_tsv <- function(x) read.delim(file.path(INPUT_DIR,x), check.names=FALSE, stringsAsFactors=FALSE)
meta <- read_tsv('latitude_sample_order.tsv')
sv <- read_tsv('latitude_variant_order.tsv')
Y <- as.matrix(read.table(file.path(INPUT_DIR,'latitude_genotype_matrix.tsv'),header=FALSE,sep='\t',check.names=FALSE))
storage.mode(Y) <- 'double'
stopifnot(nrow(Y)==nrow(meta),ncol(Y)==nrow(sv), !anyNA(Y),all(is.finite(Y)))
stopifnot(!anyDuplicated(sv$variant_key), !anyDuplicated(meta$sample_id))
if (length(unique(meta$site_id))!=11) warning('Number of sites is not 11; confirm input and sample selection')
if (nrow(Y)!=46) warning('Number of individuals differs from expected 46; verify against original LFMM data')
lat <- as.numeric(meta$lat)
if (anyNA(lat) || sd(lat)<=0) stop('Latitude missing or constant')
X <- matrix(as.numeric(scale(lat)),ncol=1,dimnames=list(meta$sample_id,'latitude'))
results <- vector('list',length(K_VALUES)); names(results) <- paste0('K',K_VALUES)
qc <- vector('list',length(K_VALUES))
for (i in seq_along(K_VALUES)) {
  k <- K_VALUES[i]
  message('Running LEA::lfmm2 for latitude with K=',k)
  model <- LEA::lfmm2(input=Y,env=X,K=k)
  # Run without automatic correction to preserve raw test z statistics
  tt <- LEA::lfmm2.test(object=model,input=Y,env=X,genomic.control=FALSE,linear=TRUE)
  z <- as.numeric(tt$zscores)
  if(length(z)!=ncol(Y)) stop('Unexpected zscore length at K=',k)
  # Same conventional median chi-square genomic-control calibration (chi-square median with 1 df)
  lambda <- median(z^2,na.rm=TRUE)/qchisq(0.5,df=1)
  # Do NOT deflate statistics if lambda<1; original workflow may differ: check original script
  lambda_used <- max(1,lambda)
  p_raw <- pchisq(z^2,df=1,lower.tail=FALSE)
  p_gc <- pchisq(z^2/lambda_used,df=1,lower.tail=FALSE)
  dat <- data.frame(variant_key=sv$variant_key,CHROM=sv$CHROM,POS=sv$POS,ID=sv$ID,
        K=k,zscore=z,p_raw=p_raw,p_genomic_control=p_gc,q_latitude_BH=p.adjust(p_gc,method='BH'),
        stringsAsFactors=FALSE)
  fn <- sprintf('05_latitude_LFMM_K%d_all_SVs.tsv',k)
  write.table(dat,file.path(OUT_DIR,fn),sep='\t',quote=FALSE,row.names=FALSE)
  results[[i]] <- dat
  qc[[i]] <- data.frame(K=k,n_individuals=nrow(Y),n_sites=length(unique(meta$site_id)),n_SVs=ncol(Y),
       genomic_inflation=lambda,lambda_used=lambda_used,significant_q05=sum(dat$q_latitude_BH<0.05,na.rm=TRUE),
       significant_q10=sum(dat$q_latitude_BH<0.1,na.rm=TRUE))
}
qc <- do.call(rbind,qc)
write.table(qc,file.path(OUT_DIR,'05_latitude_model_QC.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
merged <- results[[which(K_VALUES==PRIMARY_K)]]
for (k in setdiff(K_VALUES,PRIMARY_K)) {
 d <- results[[which(K_VALUES==k)]]
 merged[[paste0('q_K',k)]] <- d$q_latitude_BH[match(merged$variant_key,d$variant_key)]
}
merged$robust_K456_q05 <- Reduce(`&`,lapply(c(4L,5L,6L),function(k) {
 if(k==PRIMARY_K) merged$q_latitude_BH<0.05 else merged[[paste0('q_K',k)]]<0.05
}))
write.table(merged,file.path(OUT_DIR,'05_latitude_K5_with_K_sensitivity.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(merged[merged$q_latitude_BH<0.05,,drop=FALSE],file.path(OUT_DIR,'05_latitude_K5_FDR05_candidates.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
message('Complete! K5 latitude BH q<0.05: ',sum(merged$q_latitude_BH<0.05),'; robust K4/5/6: ',sum(merged$robust_K456_q05))
message('Check input filtering/imputation with the prior BIO LFMM pipeline before claiming exact comparability.')
