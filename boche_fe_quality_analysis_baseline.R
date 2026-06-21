library(haven)
library(fixest)

final <- Sys.getenv(
  "MEDICARE_FINAL_DIR",
  "E:/Data/Medicare data/AMI HF/Databases/STATA files/Final files"
)
outdir <- "r_baseline_outputs"
dir.create(outdir, showWarnings = FALSE)

d <- read_dta(file.path(final, "incident-ami-final.dta"))
demo <- read_dta(file.path(final, "other useful files", "cardiologist_demographics.dta"))

keep <- d$qualify_transfer_out == 0 &
  d$qualify_transfer_in == 0 &
  !is.na(d$npi1) & d$npi1 != "" &
  (is.na(d$ccw_ami_1y) | d$ccw_ami_1y != 1)
keep[is.na(keep)] <- FALSE
d <- d[keep, ]
d <- merge(d, demo, by = "npi1", all.x = TRUE)

d$male <- as.integer(as.numeric(d$sex) == 1)
d$race_white <- as.integer(as.numeric(d$race) == 1)
d$race_black <- as.integer(as.numeric(d$race) == 2)
d$race_asian <- as.integer(as.numeric(d$race) == 4)
d$race_hisp <- as.integer(as.numeric(d$race) == 5)

dt <- as.Date(d$qualify_admsndt)
d$admit_year_month <- as.integer(format(dt, "%Y")) * 12L + as.integer(format(dt, "%m"))
d$admit_day_week <- as.POSIXlt(dt)$wday
d$age_bins <- cut(d$qualify_age, c(64, 69, 74, 79, 84, 89, 94, Inf), labels = 1:7)
d$age_bins <- as.integer(as.character(d$age_bins))

d$card_id <- as.integer(factor(d$npi1))
d$hospital_id <- as.integer(factor(d$qualify_prvdrnum))
d$pat_vol <- ave(rep(1L, nrow(d)), d$card_id, FUN = length)
d <- d[d$pat_vol >= 5, ]

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

x <- c(
  "male", "race_white", "race_black", "race_asian", "race_hisp",
  "i(age_bins)", "i(enroll_partd)", paste0("i(", ccw, ")")
)
fml <- as.formula(paste(
  "mortality_0 ~", paste(x, collapse = " + "),
  "| card_id + hospital_id + admit_year_month + admit_day_week"
))

est <- feols(fml, data = d, notes = FALSE)
d$ehat <- resid(est)
d$yhat <- fitted(est)

fe <- fixef(est, notes = FALSE)$card_id
doc <- unique(d[c("card_id", "npi1", "pat_vol")])
doc$card_fe <- fe[match(doc$card_id, as.integer(names(fe)))]
doc <- doc[order(doc$card_id), ]

saveRDS(est, file.path(outdir, "baseline_feols_estimate.rds"))
saveRDS(d, file.path(outdir, "baseline_sample_with_residuals.rds"))
write.csv(doc, file.path(outdir, "baseline_doctor_fixed_effects.csv"), row.names = FALSE)
capture.output(summary(est), file = file.path(outdir, "baseline_summary.txt"))

cat("N =", nobs(est), "\nDoctors =", nrow(doc), "\n")
