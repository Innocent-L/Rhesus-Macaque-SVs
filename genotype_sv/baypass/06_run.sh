UTILS=/XX/XX/baypass_public/utils/baypass_utils.R

Rscript 06_make_POD.R \
    "$UTILS" \
    china46_core_mat_omega.out \
    china46_core_summary_beta_params.out \
    china46_allpoly.geno \
    04_pod \
    100000 \
    2>&1 | tee logs/04_make_POD.log
