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

drop_intercept <- function(x) {
  keep <- colnames(x) != "(Intercept)"
  x[, keep, drop = FALSE]
}

drop_first_level <- function(fml, data) {
  x <- sparse.model.matrix(fml, data = data)
  x[, -1, drop = FALSE]
}

D <- sparse.model.matrix(~ 0 + card_id, data = d)

x_fml <- as.formula(paste(
  "~ male + race_white + race_black + race_asian + race_hisp",
  "+ factor(age_bins) + factor(enroll_partd)",
  "+", paste(paste0("factor(", ccw, ")"), collapse = " + ")
))
X <- drop_intercept(sparse.model.matrix(x_fml, data = d))

H <- drop_first_level(~ 0 + hospital_id, d)
M <- drop_first_level(~ 0 + admit_year_month, d)
W <- drop_first_level(~ 0 + admit_day_week, d)

A <- cbind(D, X, H, M, W)
G <- crossprod(A)

doctor_cols <- seq_len(ncol(D))
batch <- doctor_cols[1:min(100, length(doctor_cols))]

E <- Matrix(0, nrow = ncol(A), ncol = length(batch), sparse = TRUE)

patient_weight <- as.numeric(table(d$card_id)[levels(d$card_id)])
patient_weight <- patient_weight / sum(patient_weight)

for (j in seq_along(batch)) {
  E[doctor_cols, j] <- -patient_weight
  E[batch[j], j] <- E[batch[j], j] + 1
}

Q <- tryCatch(
  solve(G, E),
  error = function(e) {
    stop(
      "Regular-case solve failed. Dropping one hospital, month, and weekday ",
      "reference was not enough to make the design full rank. Original error: ",
      conditionMessage(e),
      call. = FALSE
    )
  }
)
L <- A %*% Q
se2 <- colSums((L * L) * as.numeric(d$ehat^2))

out <- data.frame(
  card_col = colnames(A)[batch],
  card_id = sub("^card_id", "", colnames(A)[batch]),
  patient_weight = patient_weight[batch],
  centered_se = sqrt(as.numeric(se2))
)

write.csv(out, file.path(outdir, "doctor_centered_se_first_batch.csv"), row.names = FALSE)
print(summary(out$centered_se))
print(head(out))
