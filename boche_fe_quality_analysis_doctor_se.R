library(Matrix)

outdir <- "r_baseline_outputs"
d <- readRDS(file.path(outdir, "baseline_sample_with_residuals.rds"))

d <- d[!is.na(d$ehat), ]
d$card_id <- factor(d$card_id)
d$hospital_id <- factor(d$hospital_id)
d$admit_year_month <- factor(d$admit_year_month)
d$admit_day_week <- factor(d$admit_day_week)

ccw <- c(
  "ccw_alzheimer_1y", "ccw_afib_1y", "ccw_alzheimer_dementia_1y",
  "ccw_anemia_1y", "ccw_asthma_1y", "ccw_breast_cancer_1y",
  "ccw_cataract_1y", "ccw_ckd_1y", "ccw_colon_cancer_1y",
  "ccw_copd_1y", "ccw_depression_1y", "ccw_diabetes_1y",
  "ccw_endometrial_cancer_1y", "ccw_glaucoma_1y", "ccw_hip_frac_1y",
  "ccw_hyperlipidemia_1y", "ccw_hypertension_1y", "ccw_hypothyroidism_1y",
  "ccw_lung_cancer_1y", "ccw_osteoporosis_1y", "ccw_prostate_cancer_1y",
  "ccw_ra_oa_1y", "ccw_stroke_1y", "ccw_ischemic_hd_1y", "ccw_hf_1y"
)

fml <- as.formula(paste(
  "~ 0 + card_id + male + race_white + race_black + race_asian + race_hisp",
  "+ factor(age_bins) + factor(enroll_partd)",
  "+", paste(paste0("factor(", ccw, ")"), collapse = " + "),
  "+ hospital_id + admit_year_month + admit_day_week"
))

A <- sparse.model.matrix(fml, data = d)
G <- crossprod(A)

doctor_cols <- grep("^card_id", colnames(A))
batch <- doctor_cols[1:min(100, length(doctor_cols))]

E <- Matrix(0, nrow = ncol(A), ncol = length(batch), sparse = TRUE)
E[cbind(batch, seq_along(batch))] <- 1

Q <- solve(G, E)
L <- A %*% Q
se2 <- colSums((L * L) * as.numeric(d$ehat^2))

out <- data.frame(
  card_col = colnames(A)[batch],
  card_id = sub("^card_id", "", colnames(A)[batch]),
  se = sqrt(as.numeric(se2))
)

write.csv(out, file.path(outdir, "doctor_se_first_batch.csv"), row.names = FALSE)
print(summary(out$se))
print(head(out))
