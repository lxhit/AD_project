install.packages("WeightIt")
install.packages("emmeans")
install.packages("table1")
library(SparkR)
library(WeightIt)
library(emmeans)
library(table1)
sessionInfo()

AD_long <- SparkR::sql("SELECT * FROM default.bpaciotti_wan_COHORT_FINAL_CC")
print(sort(names(AD_long)))
labs <- c('BMI', 'ALP', 'ALT', 'AST', 'Protein_3020630', 'Bilirubin_Total_3024128', 'Cholesterolin_HDL_3007070', 'Cholesterol_3027114', 'Cholesterolin_LDL_3028288', 'Glucose_3004501', 'Triglyceride_3022192', 'Urea_Nitrogen_3013682', 'Calcium_3006906', 'eGFR_3049187', 'Cholesterol_Non_HDL', 'Sodium', 'Potassium', 'Creatinine', 'Heart_Rate', 'Hematocrit', 'Dia_BP_3012888', 'Sys_BP_3004249', 'Leukocytes_3000905', 'Neutrophils_3013650', 'Neutrophils_Leukocytes_3008342', 'O2_Sat_40762499', 'Platelets', 'Erythrocyte_Volume_3015182', 'MCH_3012030', 'MCV_3023599', 'MCHC_3009744', 'Hemoglobin_A1c', 'Albumin_3024561', 'C02_3015632', 'eGFR_3053283', 'Lymphocytes_Leukocytes_3037511', 'Erythrocyte_Ratio_3019897', 'Monocytes_Leukocytes_3011948', 'Eosinophils_Leukocytes_3010457')
diagnoses <- c('HCC_LIVER', 'CANCER_LIVER', 'LIVER_NAFLD', 'LIVER_NASH', 'ALCOHOLIC_RELATED_LIVER_DX', 'TOXIC_LIVER_DX', 'HEPATIC_FAILURE', 'NONINFECTIOUS_HEPATITIS', 'CIRRHOSIS', 'INFLAMMATORY_LIVER', 'ABSCESS_LIVER', 'AUTOIMMMUNE_HEPATITIS', 'LIVER_OTHER', 'LIVER_NEC', 'DIABETES', 'ALCOHOL', 'NICOTINE', 'HIV', 'CVD', 'HOT_FLASHES', 'HOT_FLASHES_MP')
AD <- SparkR::sql("SELECT * FROM default.bpaciotti_wan_PAT_SUMMARY")
AD2 <- SparkR::select(AD, AD$PERSON_ID, AD$ALZ_MIN_DT, AD$BIRTH_YEAR, AD$ADI, AD$RACE_ETHNICITY, AD$UC_SITE, AD$GENDER, AD$DEATH_DATE, AD$AGE_LAST_OBS, AD$AGE_FIRST_OBS)
AD2 <- as.data.frame(AD2)
tmp <- as.data.frame(select(AD_long, AD_long$PERSON_ID))$PERSON_ID
subjects.keep <- intersect(tmp, AD2$PERSON_ID)
AD2 <- subset(AD2, PERSON_ID %in% subjects.keep)
dim(AD2)
tmp <- sapply(labs, function(lab){
  tmp <- SparkR::select(AD_long, AD_long$PERSON_ID, AD_long[[lab]], AD_long$BIRTH_YEAR, AD_long$OBSERVATION_PERIOD_START_DATE)
  tmp <- subset(tmp, !isNull(tmp[[lab]]))
  tmp <- as.data.frame(tmp)
  tmp[[paste0("AGE_", lab)]] <- lubridate::year(tmp$OBSERVATION_PERIOD_START_DATE) - tmp$BIRTH_YEAR
  tmp <- tmp[tmp[[paste0("AGE_", lab)]] <= 69, ]
  
  # VIM for each subject
  count <- tapply(tmp[[lab]], tmp$PERSON_ID, function(x)length(x[!is.na(x)]))
  means <- tapply(tmp[[lab]], tmp$PERSON_ID, mean, na.rm = TRUE)
  sd <- tapply(tmp[[lab]], tmp$PERSON_ID, sd, na.rm = TRUE)
  means <- ifelse(count < 3, NA, means)
  sd <- ifelse(count < 3, NA, sd)
  fit <- lm(log(sd + 0.000001) ~ log(means + 0.000001), na.action = na.omit)
  x <- coef(fit)[2]
  k <- mean(means, na.rm = TRUE)^x
  VIM <- k*sd/(means^x)
  VIM.df <- data.frame(PERSON_ID = as.numeric(names(VIM)))
  VIM.df[[paste0(lab, "_VIM")]] <- VIM
  
  # Min for each subject
  min <- tapply(tmp[[lab]], tmp$PERSON_ID, min, na.rm = TRUE)
  min.df <- data.frame(PERSON_ID = as.numeric(names(min)))
  min.df[[paste0(lab, "_MIN")]] <- min
  
  # Max for each subject
  max <- tapply(tmp[[lab]], tmp$PERSON_ID, max, na.rm = TRUE)
  max.df <- data.frame(PERSON_ID = as.numeric(names(max)))
  max.df[[paste0(lab, "_MAX")]] <- max
  head(max.df)
  
  # latest observation for each subject
  tmp2 <- dplyr::arrange(tmp, PERSON_ID, desc(OBSERVATION_PERIOD_START_DATE))
  tmp2 <- dplyr::group_by(tmp2, PERSON_ID)
  tmp2 <- dplyr::slice(tmp2, 1)
  tmp2 <- tmp2[,c(lab, "PERSON_ID", paste0("AGE_", lab))]
  names(tmp2)[c(1, 3)] <- c(paste0(lab, "_LATEST"), paste0("AGE_LATEST_", lab))
  head(tmp2)
  
  AD2 <<- dplyr::left_join(AD2, tmp2, by = "PERSON_ID")
  AD2 <<- dplyr::left_join(AD2, min.df, by = "PERSON_ID")
  AD2 <<- dplyr::left_join(AD2, max.df, by = "PERSON_ID")
  AD2 <<- dplyr::left_join(AD2, VIM.df, by = "PERSON_ID")
  return(NULL)
})
names(AD2)
tmp <- sapply(diagnoses, function(dx){
  tmp <- select(AD_long, AD_long$PERSON_ID, AD_long[[dx]], AD_long$BIRTH_YEAR, AD_long$OBSERVATION_PERIOD_START_DATE)
  tmp <- subset(tmp, tmp[[dx]] == 1)
  tmp <- as.data.frame(tmp)
  tmp[[paste0("AGE_", dx)]] <- lubridate::year(tmp$OBSERVATION_PERIOD_START_DATE) - tmp$BIRTH_YEAR
  tmp <- tmp[tmp[[paste0("AGE_", dx)]] <= 69, ]
  # Take the latest observation for each subject
  tmp2 <- dplyr::arrange(tmp, PERSON_ID, desc(OBSERVATION_PERIOD_START_DATE))
  tmp2 <- dplyr::group_by(tmp2, PERSON_ID)
  tmp2 <- dplyr::slice(tmp2, 1)
  tmp2 <- tmp2[,c(dx, "PERSON_ID")]
  AD2 <<- dplyr::left_join(AD2, tmp2, by = "PERSON_ID")
  return(NULL)
})
%md
## Excluded anyone who died or had AD prior to age 70
AD2$AD_DX_AGE <- lubridate::year(AD2$ALZ_MIN_DT) - AD2$BIRTH_YEAR
AD2$AD_DX_IND <- ifelse(is.na(AD2$AD_DX_AGE), 0, 1)
AD2$DEATH_AGE <- lubridate::year(AD2$DEATH_DATE) - AD2$BIRTH_YEAR
AD2$DEATH_IND <- ifelse(is.na(AD2$DEATH_AGE), 0, 1)
bad <- which((!is.na(AD2$AD_DX_AGE) & AD2$AD_DX_AGE <= 70) | (!is.na(AD2$DEATH_AGE) & AD2$DEATH_AGE <= 70))
AD2 <- AD2[-bad,]
print(dim(AD2))
%md
## Exclude anyone whose age at first visit is 69 or older
bad <- which(AD2$AGE_FIRST_OBS >= 69)
AD2 <- AD2[-bad,]
dim(AD2)
%md
### Min and max variables with normal ranges split by those, other min/max variables split into tertiles, VIMs split into quintiles
# Categorize variables
# eGFR normal range >= 60
AD2$eGFR_3049187_MIN_cat <- cut(AD2$eGFR_3049187_MIN, breaks = c(-Inf, 60, Inf), right = FALSE)
levels(AD2$eGFR_3049187_MIN_cat) <- c("Low", "Normal")
AD2$eGFR_3053283_MIN_cat <- cut(AD2$eGFR_3053283_MIN, breaks = c(-Inf, 60, Inf), right = FALSE)
levels(AD2$eGFR_3053283_MIN_cat) <- c("Low", "Normal")
AD2$eGFR_MIN_cat <- NA
tmp <- which(AD2$eGFR_3049187_MIN_cat == "Low" | AD2$eGFR_3053283_MIN_cat == "Low")
AD2$eGFR_MIN_cat[tmp] <- "Low"
tmp <- which(AD2$eGFR_3049187_MIN_cat == "Normal" & AD2$eGFR_3053283_MIN_cat == "Normal")
AD2$eGFR_MIN_cat[tmp] <- "Normal"
tmp <- which((AD2$eGFR_3049187_MIN_cat == "Normal" & is.na(AD2$eGFR_3053283_MIN_cat)) | (is.na(AD2$eGFR_3049187_MIN_cat) & AD2$eGFR_3053283_MIN_cat == "Normal"))
AD2$eGFR_MIN_cat[tmp] <- "Normal"
# Erythrocyte_Ratio_3019897 normal range is 0 to 15 mm/hr in men. 0 to 20 mm/hr in women
male.cat <- cut(AD2$Erythrocyte_Ratio_3019897_MIN, breaks = c(-Inf, 15, Inf))
levels(male.cat) <- c("Normal", "High")
female.cat <- cut(AD2$Erythrocyte_Ratio_3019897_MIN, breaks = c(-Inf, 20, Inf))
levels(female.cat) <- c("Normal", "High")
AD2$Erythrocyte_Ratio_3019897_MIN_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$Erythrocyte_Ratio_3019897_MIN_cat) <- c("Normal", "High")
# Sodium normal range is 136 - 145 mmol/L
AD2$Sodium_MIN_cat <- cut(AD2$Sodium_MIN, breaks = c(-Inf, 135.9, 145, Inf))
levels(AD2$Sodium_MIN_cat) <- c("Low", "Normal", "High")
# Urea_Nitrogen_3013682 normal range is 6 - 20 mg/dL
AD2$Urea_Nitrogen_3013682_MIN_cat <- cut(AD2$Urea_Nitrogen_3013682_MIN, breaks = c(-Inf, 5.99, 20, Inf))
levels(AD2$Urea_Nitrogen_3013682_MIN_cat) <- c("Low", "Normal", "High")
# Heart_Rate normal range is 60 to 100
AD2$Heart_Rate_MIN_cat <- cut(AD2$Heart_Rate_MIN, breaks = c(-Inf, 59.99, 100, Inf))
levels(AD2$Heart_Rate_MIN_cat) <- c("Low", "Normal", "High")
# Neutrophils_Leukocytes_3008342
# AST normal range is 9 - 32 U/L for female 10 - 40 U/L for male
male.cat <- cut(AD2$AST_MIN, breaks = c(-Inf, 9.99, 40, Inf))
levels(male.cat) <- c("Low", "Normal", "High")
female.cat <- cut(AD2$AST_MIN, breaks = c(-Inf, 8.99, 32, Inf))
levels(female.cat) <- c("Low", "Normal", "High")
AD2$AST_MIN_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$AST_MIN_cat) <- c("Low", "Normal", "High")
# Albumin_3024561 normal range is 4.0 - 4.9 g/dL
AD2$Albumin_3024561_MIN_cat <- cut(AD2$Albumin_3024561_MIN, breaks = c(-Inf, 3.99, 4.9, Inf))
levels(AD2$Albumin_3024561_MIN_cat) <- c("Low", "Normal", "High")
# ALP normal range is 44-147 U/L
AD2$ALP_MIN_cat <- cut(AD2$ALP_MIN, breaks = c(-Inf, 43.99, 147, Inf))
levels(AD2$ALP_MIN_cat) <- c("Low", "Normal", "High")
# Monocytes_Leukocytes_3011948
# Glucose_3004501 normal range is less than 140 mg/dL
AD2$Glucose_3004501_MIN_cat <- cut(AD2$Glucose_3004501_MIN, breaks = c(-Inf, 140, Inf), right = FALSE)
levels(AD2$Glucose_3004501_MIN_cat) <- c("Normal", "High")
# BP
AD2$BP_MIN_cat <- ifelse(AD2$Sys_BP_3004249_MIN < 90 & AD2$Dia_BP_3012888_MIN < 60,
                         "Low",
                         ifelse(AD2$Sys_BP_3004249_MIN < 120 & AD2$Dia_BP_3012888_MIN < 80, 
                                "Normal",
                                ifelse(AD2$Sys_BP_3004249_MIN < 140 & AD2$Dia_BP_3012888_MIN < 90, 
                                       "HTN I",
                                       ifelse(AD2$Sys_BP_3004249_MIN < 180 & AD2$Dia_BP_3012888_MIN < 120,
                                              "HTN II",
                                              ifelse(!is.na(AD2$Sys_BP_3004249_MIN) & !is.na(AD2$Dia_BP_3012888_MIN), 
                                                     "Hypertensive Crisis", NA)))))
AD2$BP_MIN_cat <- factor(AD2$BP_MIN_cat, levels = c("Low", "Normal", "HTN I", "HTN II", "Hypertensive Crisis"))
# Lymphocytes_Leukocytes_3037511
# Triglyceride_3022192 normal range is Less than 150mg/dL
AD2$Triglyceride_3022192_MIN_cat <- cut(AD2$Triglyceride_3022192_MIN, breaks = c(-Inf, 150, Inf), right = FALSE)
levels(AD2$Triglyceride_3022192_MIN_cat) <- c("Normal", "High")
# Neutrophils_3013650 normal range is 1.8 - 7.7 K/MM3
AD2$Neutrophils_3013650_MIN_cat <- cut(AD2$Neutrophils_3013650_MIN, breaks = c(-Inf, 1.79, 7.7, Inf))
levels(AD2$Neutrophils_3013650_MIN_cat) <- c("Low", "Normal", "High")
# Platelets normal range is 130 - 400 K/MM3
AD2$Platelets_MIN_cat <- cut(AD2$Platelets_MIN, breaks = c(-Inf, 129.9, 400, Inf))
levels(AD2$Platelets_MIN_cat) <- c("Low", "Normal", "High")
# Hematocrit normal range is 41-50 for men and 36-48 for women
male.cat <- cut(AD2$Hematocrit_MIN, breaks = c(-Inf, 40.99, 50, Inf), right = TRUE)
levels(male.cat) <- c("Low", "Normal", "High")
female.cat <- cut(AD2$Hematocrit_MIN, breaks = c(-Inf, 35.99, 48, Inf), right = TRUE)
levels(female.cat) <- c("Low", "Normal", "High")
AD2$Hematocrit_MIN_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$Hematocrit_MIN_cat) = c("Low", "Normal", "High")
# Erythrocyte_Volume_3015182
# Cholesterolin_HDL_3007070 normal range is 40mg/dL or higher for women and 50mg/dL or higher for men
male.cat <- cut(AD2$Cholesterolin_HDL_3007070_MIN, breaks = c(-Inf, 50, Inf), right = FALSE)
levels(male.cat) <- c("Low", "Normal", "High")
female.cat <- cut(AD2$Cholesterolin_HDL_3007070_MIN, breaks = c(-Inf, 40, Inf), right = FALSE)
levels(male.cat) <- c("Low", "Normal", "High")
AD2$Cholesterolin_HDL_3007070_MIN_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$Cholesterolin_HDL_3007070_MIN_cat) = c("Low", "Normal")
# Leukocytes_3000905
# c02_3015632 normal range is 22-29 mmol/L
AD2$C02_3015632_MIN_cat <- cut(AD2$C02_3015632_MIN, breaks = c(-Inf, 21.99, 29, Inf))
levels(AD2$C02_3015632_MIN_cat) <- c("Low", "Normal", "High")
# Hemoglobin_A1c Nomral: <= 5.7 Prediabetes: > 5.7 & < 6.4 Diabetes: >= 6.4
AD2$Hemoglobin_A1c_MIN_cat <- cut(AD2$Hemoglobin_A1c_MIN, breaks = c(-Inf, 5.7, 6.39, Inf))
levels(AD2$Hemoglobin_A1c_MIN_cat) <- c("Normal", "Prediabetes", "Diabetes")
# Eosinophils_Leukocytes_3010457
# Bilirubin_Total_3024128 normal range is <= 1.2 mg/dL
AD2$Bilirubin_Total_3024128_MIN_cat <- cut(AD2$Bilirubin_Total_3024128_MIN, breaks = c(-Inf, 1.2, Inf), right = FALSE)
levels(AD2$Bilirubin_Total_3024128_MIN_cat) <- c("Normal", "High")
# Protein_3020630 normal range is 6.6 - 8.7 g/dL
AD2$Protein_3020630_MIN_cat <- cut(AD2$Protein_3020630_MIN, breaks = c(-Inf, 6.59, 8.7, Inf))
levels(AD2$Protein_3020630_MIN_cat) <- c("Low", "Normal", "High")
# ALT normal range is 4 to 36 U/L
AD2$ALT_MIN_cat <- cut(AD2$ALT_MIN, breaks = c(-Inf, 3.99, 36, Inf))
levels(AD2$ALT_MIN_cat) <- c("Low", "Normal", "High")
# Calcium_300690 normal range is 8.6 - 10.0 mg/dL
AD2$Calcium_3006906_MIN_cat <- cut(AD2$Calcium_3006906_MIN, breaks = c(-Inf, 8.59, 10, Inf))
levels(AD2$Calcium_3006906_MIN_cat) <- c("Low", "Normal", "High")
# Normal range for Cholesterol_3027114 is 125 to 200mg/dL
AD2$Cholesterol_3027114_MIN_cat <- cut(AD2$Cholesterol_3027114_MIN, breaks = c(-Inf, 125, 200, Inf))
levels(AD2$Cholesterol_3027114_MIN_cat) <- c("Low", "Normal", "High")
# Normal range for Cholesterolin_LDL_3028288 is < 100
AD2$Cholesterolin_LDL_3028288_MIN_cat <- cut(AD2$Cholesterolin_LDL_3028288_MIN, breaks = c(-Inf, 100, Inf))
levels(AD2$Cholesterolin_LDL_3028288_MIN_cat) <- c("Normal", "High")
# Normal range for Cholesterol_Non_HDL is less than 130
AD2$Cholesterol_Non_HDL_MIN_cat <- cut(AD2$Cholesterol_Non_HDL_MIN, breaks = c(-Inf, 130, Inf))
levels(AD2$Cholesterol_Non_HDL_MIN_cat) <- c("Normal", "High")
# Normal range for potassium is 3.4 - 5.1 mmol/L
AD2$Potassium_MIN_cat <- cut(AD2$Potassium_MIN, breaks = c(-Inf, 3.4, 5.1, Inf))
levels(AD2$Potassium_MIN_cat) <- c("Low", "Normal", "High")
# Normal range for creatinine is Female	0.59 - 1.04 mg/dL Male	0.74 - 1.35 mg/dL
male.cat <- cut(AD2$Creatinine_MIN, breaks = c(-Inf, 0.74, 1.35, Inf), right = TRUE)
levels(male.cat) <- c("Low", "Normal", "High")
female.cat <- cut(AD2$Creatinine_MIN, breaks = c(-Inf, 0.59, 1.04, Inf), right = TRUE)
levels(female.cat) <- c("Low", "Normal", "High")
AD2$Creatinine_MIN_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$Creatinine_MIN_cat) = c("Low", "Normal", "High")
# Everything else split by tertiles
other.vars <-  c("Leukocytes_3000905_MIN", "Neutrophils_Leukocytes_3008342_MIN", "O2_Sat_40762499_MIN", "Erythrocyte_Volume_3015182_MIN", "MCH_3012030_MIN", "MCV_3023599_MIN", "MCHC_3009744_MIN", "Lymphocytes_Leukocytes_3037511_MIN", "Monocytes_Leukocytes_3011948_MIN", "Eosinophils_Leukocytes_3010457_MIN")
for (var in other.vars){
  newvar <- paste0(var, "_cat")
  AD2[[newvar]] <- cut(AD2[[var]], breaks = c(-Inf, quantile(AD2[[var]], c(1/3, 2/3), na.rm = TRUE), Inf))
  levels(AD2[[newvar]]) <- c("Low", "Normal", "High")
}
names(AD2)
# Cleaning BMI
AD2$BMI_MIN <- ifelse(AD2$BMI_MIN < 12, NA, AD2$BMI_MIN)
# BMI_MIN_MIN_cat_CDC
AD2$BMI_MIN_cat_CDC <- cut(AD2$BMI_MIN, breaks = c(-Inf, 18.5, 25, 30, Inf), right = FALSE)
levels(AD2$BMI_MIN_cat_CDC) <- c("Underweight", "Normal", "Overweight", "Obese")
# BMI_MIN_cat_RACE
tmp <-  cut(AD2$BMI_MIN, breaks = c(-Inf, 18.5, 23, 27.5, Inf), right = FALSE)
levels(tmp) <- c("Underweight", "Normal", "Overweight", "Obese")
AD2$BMI_MIN_cat_RACE <- factor(ifelse(AD2$RACE_ETHNICITY == "Asian", as.character(tmp), as.character(AD2$BMI_MIN_cat_CDC)), levels = c("Underweight", "Normal", "Overweight", "Obese"))
# Categorize variables
# eGFR normal range >= 60
AD2$eGFR_3049187_MAX_cat <- cut(AD2$eGFR_3049187_MAX, breaks = c(-Inf, 60, Inf), right = FALSE)
levels(AD2$eGFR_3049187_MAX_cat) <- c("Low", "Normal")
AD2$eGFR_3053283_MAX_cat <- cut(AD2$eGFR_3053283_MAX, breaks = c(-Inf, 60, Inf), right = FALSE)
levels(AD2$eGFR_3053283_MAX_cat) <- c("Low", "Normal")
AD2$eGFR_MAX_cat <- NA
tmp <- which(AD2$eGFR_3049187_MAX_cat == "Low" | AD2$eGFR_3053283_MAX_cat == "Low")
AD2$eGFR_MAX_cat[tmp] <- "Low"
tmp <- which(AD2$eGFR_3049187_MAX_cat == "Normal" & AD2$eGFR_3053283_MAX_cat == "Normal")
AD2$eGFR_MAX_cat[tmp] <- "Normal"
tmp <- which((AD2$eGFR_3049187_MAX_cat == "Normal" & is.na(AD2$eGFR_3053283_MAX_cat)) | (is.na(AD2$eGFR_3049187_MAX_cat) & AD2$eGFR_3053283_MAX_cat == "Normal"))
AD2$eGFR_MAX_cat[tmp] <- "Normal"
# Erythrocyte_Ratio_3019897 normal range is 0 to 15 mm/hr in men. 0 to 20 mm/hr in women
male.cat <- cut(AD2$Erythrocyte_Ratio_3019897_MAX, breaks = c(-Inf, 15, Inf))
levels(male.cat) <- c("Normal", "High")
female.cat <- cut(AD2$Erythrocyte_Ratio_3019897_MAX, breaks = c(-Inf, 20, Inf))
levels(female.cat) <- c("Normal", "High")
AD2$Erythrocyte_Ratio_3019897_MAX_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$Erythrocyte_Ratio_3019897_MAX_cat) <- c("Normal", "High")
# Sodium normal range is 136 - 145 mmol/L
AD2$Sodium_MAX_cat <- cut(AD2$Sodium_MAX, breaks = c(-Inf, 135.9, 145, Inf))
levels(AD2$Sodium_MAX_cat) <- c("Low", "Normal", "High")
# Urea_Nitrogen_3013682 normal range is 6 - 20 mg/dL
AD2$Urea_Nitrogen_3013682_MAX_cat <- cut(AD2$Urea_Nitrogen_3013682_MAX, breaks = c(-Inf, 5.99, 20, Inf))
levels(AD2$Urea_Nitrogen_3013682_MAX_cat) <- c("Low", "Normal", "High")
# Heart_Rate normal range is 60 to 100
AD2$Heart_Rate_MAX_cat <- cut(AD2$Heart_Rate_MAX, breaks = c(-Inf, 59.99, 100, Inf))
levels(AD2$Heart_Rate_MAX_cat) <- c("Low", "Normal", "High")
# Neutrophils_Leukocytes_3008342
# AST normal range is 9 - 32 U/L for female 10 - 40 U/L for male
male.cat <- cut(AD2$AST_MAX, breaks = c(-Inf, 9.99, 40, Inf))
levels(male.cat) <- c("Low", "Normal", "High")
female.cat <- cut(AD2$AST_MAX, breaks = c(-Inf, 8.99, 32, Inf))
levels(female.cat) <- c("Low", "Normal", "High")
AD2$AST_MAX_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$AST_MAX_cat) <- c("Low", "Normal", "High")
# Albumin_3024561 normal range is 4.0 - 4.9 g/dL
AD2$Albumin_3024561_MAX_cat <- cut(AD2$Albumin_3024561_MAX, breaks = c(-Inf, 3.99, 4.9, Inf))
levels(AD2$Albumin_3024561_MAX_cat) <- c("Low", "Normal", "High")
# ALP normal range is 44-147 U/L
AD2$ALP_MAX_cat <- cut(AD2$ALP_MAX, breaks = c(-Inf, 43.99, 147, Inf))
levels(AD2$ALP_MAX_cat) <- c("Low", "Normal", "High")
# Monocytes_Leukocytes_3011948
# Glucose_3004501 normal range is less than 140 mg/dL
AD2$Glucose_3004501_MAX_cat <- cut(AD2$Glucose_3004501_MAX, breaks = c(-Inf, 140, Inf), right = FALSE)
levels(AD2$Glucose_3004501_MAX_cat) <- c("Normal", "High")
# BP
AD2$BP_MAX_cat <- ifelse(AD2$Sys_BP_3004249_MAX < 90 & AD2$Dia_BP_3012888_MAX < 60,
                         "Low",
                         ifelse(AD2$Sys_BP_3004249_MAX < 120 & AD2$Dia_BP_3012888_MAX < 80, 
                                "Normal",
                                ifelse(AD2$Sys_BP_3004249_MAX < 140 & AD2$Dia_BP_3012888_MAX < 90, 
                                       "HTN I",
                                       ifelse(AD2$Sys_BP_3004249_MAX < 180 & AD2$Dia_BP_3012888_MAX < 120,
                                              "HTN II",
                                              ifelse(!is.na(AD2$Sys_BP_3004249_MAX) & !is.na(AD2$Dia_BP_3012888_MAX), 
                                                     "Hypertensive Crisis", NA)))))
AD2$BP_MAX_cat <- factor(AD2$BP_MAX_cat, levels = c("Low", "Normal", "HTN I", "HTN II", "Hypertensive Crisis"))
# Lymphocytes_Leukocytes_3037511
# Triglyceride_3022192 normal range is Less than 150mg/dL
AD2$Triglyceride_3022192_MAX_cat <- cut(AD2$Triglyceride_3022192_MAX, breaks = c(-Inf, 150, Inf), right = FALSE)
levels(AD2$Triglyceride_3022192_MAX_cat) <- c("Normal", "High")
# Neutrophils_3013650 normal range is 1.8 - 7.7 K/MM3
AD2$Neutrophils_3013650_MAX_cat <- cut(AD2$Neutrophils_3013650_MAX, breaks = c(-Inf, 1.79, 7.7, Inf))
levels(AD2$Neutrophils_3013650_MAX_cat) <- c("Low", "Normal", "High")
# Platelets normal range is 130 - 400 K/MM3
AD2$Platelets_MAX_cat <- cut(AD2$Platelets_MAX, breaks = c(-Inf, 129.9, 400, Inf))
levels(AD2$Platelets_MAX_cat) <- c("Low", "Normal", "High")
# Hematocrit normal range is 41-50 for men and 36-48 for women
male.cat <- cut(AD2$Hematocrit_MAX, breaks = c(-Inf, 40.99, 50, Inf), right = TRUE)
levels(male.cat) <- c("Low", "Normal", "High")
female.cat <- cut(AD2$Hematocrit_MAX, breaks = c(-Inf, 35.99, 48, Inf), right = TRUE)
levels(female.cat) <- c("Low", "Normal", "High")
AD2$Hematocrit_MAX_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$Hematocrit_MAX_cat) = c("Low", "Normal", "High")
# Erythrocyte_Volume_3015182
# Cholesterolin_HDL_3007070 normal range is 40mg/dL or higher for women and 50mg/dL or higher for men
male.cat <- cut(AD2$Cholesterolin_HDL_3007070_MAX, breaks = c(-Inf, 50, Inf), right = FALSE)
levels(male.cat) <- c("Low", "Normal", "High")
female.cat <- cut(AD2$Cholesterolin_HDL_3007070_MAX, breaks = c(-Inf, 40, Inf), right = FALSE)
levels(male.cat) <- c("Low", "Normal", "High")
AD2$Cholesterolin_HDL_3007070_MAX_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$Cholesterolin_HDL_3007070_MAX_cat) = c("Low", "Normal")
# Leukocytes_3000905
# c02_3015632 normal range is 22-29 mmol/L
AD2$C02_3015632_MAX_cat <- cut(AD2$C02_3015632_MAX, breaks = c(-Inf, 21.99, 29, Inf))
levels(AD2$C02_3015632_MAX_cat) <- c("Low", "Normal", "High")
# Hemoglobin_A1c Nomral: <= 5.7 Prediabetes: > 5.7 & < 6.4 Diabetes: >= 6.4
AD2$Hemoglobin_A1c_MAX_cat <- cut(AD2$Hemoglobin_A1c_MAX, breaks = c(-Inf, 5.7, 6.39, Inf))
levels(AD2$Hemoglobin_A1c_MAX_cat) <- c("Normal", "Prediabetes", "Diabetes")
# Eosinophils_Leukocytes_3010457
# Bilirubin_Total_3024128 normal range is <= 1.2 mg/dL
AD2$Bilirubin_Total_3024128_MAX_cat <- cut(AD2$Bilirubin_Total_3024128_MAX, breaks = c(-Inf, 1.2, Inf), right = FALSE)
levels(AD2$Bilirubin_Total_3024128_MAX_cat) <- c("Normal", "High")
# Protein_3020630 normal range is 6.6 - 8.7 g/dL
AD2$Protein_3020630_MAX_cat <- cut(AD2$Protein_3020630_MAX, breaks = c(-Inf, 6.59, 8.7, Inf))
levels(AD2$Protein_3020630_MAX_cat) <- c("Low", "Normal", "High")
# ALT normal range is 4 to 36 U/L
AD2$ALT_MAX_cat <- cut(AD2$ALT_MAX, breaks = c(-Inf, 3.99, 36, Inf))
levels(AD2$ALT_MAX_cat) <- c("Low", "Normal", "High")
# Calcium_300690 normal range is 8.6 - 10.0 mg/dL
AD2$Calcium_3006906_MAX_cat <- cut(AD2$Calcium_3006906_MAX, breaks = c(-Inf, 8.59, 10, Inf))
levels(AD2$Calcium_3006906_MAX_cat) <- c("Low", "Normal", "High")
# Normal range for Cholesterol_3027114 is 125 to 200mg/dL
AD2$Cholesterol_3027114_MAX_cat <- cut(AD2$Cholesterol_3027114_MAX, breaks = c(-Inf, 125, 200, Inf))
levels(AD2$Cholesterol_3027114_MAX_cat) <- c("Low", "Normal", "High")
# Normal range for Cholesterolin_LDL_3028288 is < 100
AD2$Cholesterolin_LDL_3028288_MAX_cat <- cut(AD2$Cholesterolin_LDL_3028288_MAX, breaks = c(-Inf, 100, Inf))
levels(AD2$Cholesterolin_LDL_3028288_MAX_cat) <- c("Normal", "High")
# Normal range for Cholesterol_Non_HDL is less than 130
AD2$Cholesterol_Non_HDL_MAX_cat <- cut(AD2$Cholesterol_Non_HDL_MAX, breaks = c(-Inf, 130, Inf))
levels(AD2$Cholesterol_Non_HDL_MAX_cat) <- c("Normal", "High")
# Normal range for potassium is 3.4 - 5.1 mmol/L
AD2$Potassium_MAX_cat <- cut(AD2$Potassium_MAX, breaks = c(-Inf, 3.4, 5.1, Inf))
levels(AD2$Potassium_MAX_cat) <- c("Low", "Normal", "High")
# Normal range for creatinine is Female	0.59 - 1.04 mg/dL Male	0.74 - 1.35 mg/dL
male.cat <- cut(AD2$Creatinine_MAX, breaks = c(-Inf, 0.74, 1.35, Inf), right = TRUE)
levels(male.cat) <- c("Low", "Normal", "High")
female.cat <- cut(AD2$Creatinine_MAX, breaks = c(-Inf, 0.59, 1.04, Inf), right = TRUE)
levels(female.cat) <- c("Low", "Normal", "High")
AD2$Creatinine_MAX_cat <- factor(ifelse(AD2$GENDER == "FEMALE", female.cat, male.cat))
levels(AD2$Creatinine_MAX_cat) = c("Low", "Normal", "High")
# Everything else split by tertiles
other.vars <-  c("Leukocytes_3000905_MAX", "Neutrophils_Leukocytes_3008342_MAX", "O2_Sat_40762499_MAX", "Erythrocyte_Volume_3015182_MAX", "MCH_3012030_MAX", "MCV_3023599_MAX", "MCHC_3009744_MAX", "Lymphocytes_Leukocytes_3037511_MAX", "Monocytes_Leukocytes_3011948_MAX", "Eosinophils_Leukocytes_3010457_MAX")
for (var in other.vars){
  newvar <- paste0(var, "_cat")
  AD2[[newvar]] <- cut(AD2[[var]], breaks = c(-Inf, quantile(AD2[[var]], c(1/3, 2/3), na.rm = TRUE), Inf))
  levels(AD2[[newvar]]) <- c("Low", "Normal", "High")
}
names(AD2)
# Cleaning BMI
AD2$BMI_MAX <- ifelse(AD2$BMI_MAX < 12, NA, AD2$BMI_MAX)
# BMI_MAX_MAX_cat_CDC
AD2$BMI_MAX_cat_CDC <- cut(AD2$BMI_MAX, breaks = c(-Inf, 18.5, 25, 30, Inf), right = FALSE)
levels(AD2$BMI_MAX_cat_CDC) <- c("Underweight", "Normal", "Overweight", "Obese")
# BMI_MAX_cat_RACE
tmp <-  cut(AD2$BMI_MAX, breaks = c(-Inf, 18.5, 23, 27.5, Inf), right = FALSE)
levels(tmp) <- c("Underweight", "Normal", "Overweight", "Obese")
AD2$BMI_MAX_cat_RACE <- factor(ifelse(AD2$RACE_ETHNICITY == "Asian", as.character(tmp), as.character(AD2$BMI_MAX_cat_CDC)), levels = c("Underweight", "Normal", "Overweight", "Obese"))
# Categorize _VIM variables
vimvars <- grep("VIM", names(AD2), value = TRUE)
for (var in vimvars){
  newvar <- paste0(var, "_cat")
  AD2[[newvar]] <- cut(AD2[[var]], breaks = c(-Inf, quantile(AD2[[var]], c(1/5, 2/5, 3/5, 4/5), na.rm = TRUE), Inf))
  levels(AD2[[newvar]]) <- c("q1", "q2", "q3", "q4", "q5")
}
AD2 <- subset(AD2, !is.na(ADI) & !(GENDER %in% c("OTHER", "UNKNOWN")) & !(RACE_ETHNICITY %in% c("Unknown", "Other Race")))
print(dim(AD2))
%md
## Cause specific survival analysis of time to AD 
Death treated as censoring event
library(survival)
# Derive first of AD, death, lost to followup
AD2$ftime <- apply(dplyr::select(AD2, AGE_LAST_OBS, AD_DX_AGE, DEATH_AGE), 1, min, na.rm = TRUE) 
# For patients for whom death is last observation, need to make sure they are flagged as dead
# For patients for whom AD dx is last observation, need to flag that
AD2$fstatus <- apply(dplyr::select(AD2, AGE_LAST_OBS, AD_DX_AGE, DEATH_AGE), 1, function(x){
  min <- min(x, na.rm = TRUE)
  return(max(which(x == min)))
})
AD2$fstatus <- AD2$fstatus - 1
print(table(AD2$fstatus))
%md
## Table 1


tabdat <- AD2
vars.use <- c("Age", "Gender", "Race", "ADI")
tabdat$Status <- factor(tabdat$fstatus, levels = c(1, 2, 0))
levels(tabdat$Status) <- c("AD", "Died Without AD", "Lost to Followup")
tabdat$Race <- ifelse(tabdat$RACE_ETHNICITY %in% c("Multirace", "American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander"), "Other", tabdat$RACE_ETHNICITY)
tabdat$Race <- forcats::fct_relevel(tabdat$Race, "Other", after = Inf)
tabdat$Age <- tabdat$AGE_LAST_OBS
tabdat$Gender <- factor(tabdat$GENDER)
levels(tabdat$Gender) <- c("Female", "Male")
form <- paste0("~", paste(vars.use, collapse = "+"), "|Status")
tab <- as.data.frame(table1(as.formula(form), data = tabdat))
display(tab)
df_tab <- createDataFrame(tab)
createOrReplaceTempView(df_tab, "df_tab")
%md
## Diagnoses
Models were adjusted for gender, race, UC site, and ADI
Models are NOT weighted for missingness since there is no distinction in the data between missing and not having a diagnosis.
dx.results <- lapply(diagnoses, function(var){
  df <- AD2
  df$x <- !is.na(df[[var]])
  df$x <- ifelse(grepl("HOT_FLASHES", var) & df$GENDER == "Male", NA, df$x) # set hot flashes to missing for men
  df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
  mod <- coxph(Surv(ftime, cenvar) ~x + GENDER + RACE_ETHNICITY + UC_SITE + ADI, data = df) 
  coefs <- summary(mod)$coef
  out <- as.data.frame(coefs[grep("x", rownames(coefs)), c(2, 3, 5), drop = FALSE])
  out$Variable <- paste0(var, ": Yes vs. No")
  out <- dplyr::select(out, Variable, everything())
  return(out)
})
dx.results2 <- as.data.frame(do.call("rbind", dx.results))
dx.results2$adj.P.Val <- p.adjust(dx.results2$P, "bonferroni")
dx.results2$Sig <- ifelse(dx.results2$adj.P.Val < 0.05, "Significant", "")
display(dx.results2)
df_tab <- createDataFrame(dx.results2)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/dx_cause_specific.csv")
%md
## Effect of alcohol on AD by race
Excluded AI/AK native and multirace subjects since numbers too small 
df <- subset(AD2, !(RACE_ETHNICITY %in% c("American Indian or Alaska Native", "Multirace")))
var <- "ALCOHOL"
df$x <- !is.na(df[[var]])
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~x*(RACE_ETHNICITY + GENDER) + UC_SITE + ADI, data = df) 
pairs0 <- pairs(emmeans(mod, ~x|RACE_ETHNICITY), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, RACE_ETHNICITY, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = "RACE_ETHNICITY")
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp$contrast <- "ALCOHOL Yes - No"
tmp <- dplyr::select(tmp, RACE_ETHNICITY, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/alcohol_by_race.csv")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp$contrast <- "Alcohol Yes - No"
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("White", contrast1))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/alcohol_race_interaction.csv")
%md 
## Effects of alcohol on AD by gender
pairs0 <- pairs(emmeans(mod, ~x|GENDER), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, GENDER, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = "GENDER")
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp$contrast <- "ALCOHOL Yes - No"
tmp <- dplyr::select(tmp, GENDER, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/alcohol_by_gender.csv")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp$contrast <- "Alcohol Yes - No"
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/alcohol_gender_interaction.csv")
%md
## Effects of NONINFECTIOUS_HEPATITIS on AD by race/gender could not be analyzed because numbers are too small
df <- AD2
var <- "NONINFECTIOUS_HEPATITIS"
df$x <- !is.na(df[[var]])
print(table(df$x, df$RACE_ETHNICITY))
print(table(df$x, df$GENDER))
%md
## Labs
Models were adjusted for gender, race, UC site, and ADI
# Clean up some small categories
AD2$Sodium_MIN_cat <- ifelse(AD2$Sodium_MIN_cat == "High", NA, as.character(AD2$Sodium_MIN_cat))
AD2$Protein_3020630_MIN_cat <- ifelse(AD2$Protein_3020630_MIN_cat == "High", NA, as.character(AD2$Protein_3020630_MIN_cat))
AD2$Urea_Nitrogen_3013682_MAX_cat <- ifelse(AD2$Urea_Nitrogen_3013682_MAX_cat == "Low", NA, as.character(AD2$Urea_Nitrogen_3013682_MAX_cat))
AD2$AST_MAX_cat <- ifelse(AD2$AST_MAX_cat == "Low", NA, as.character(AD2$AST_MAX_cat))

catvars <- grep("_MIN_cat|_MAX_cat", names(AD2), value = TRUE)
catvars <- setdiff(catvars, c("BMI_MIN_cat_CDC", "BMI_MIN_cat_RACE", "BMI_MAX_cat_CDC", "BMI_MAX_cat_RACE", "eGFR_3049187_MIN_cat", "eGFR_3053283_MIN_cat", "eGFR_3049187_MAX_cat", "eGFR_3053283_MAX_cat"))
lab.results <- lapply(catvars, function(var){
  df <- AD2
  df$x <- df[[var]]
  # Model missingness
  df$nonmissing <- !is.na(df$x)
  weights <- weightit(nonmissing ~ RACE_ETHNICITY + GENDER + AGE_LAST_OBS + ADI + UC_SITE, data = df, estimand = "ATE", method = "ps")
  df$w <- weights$weights
  df <- df[df$nonmissing,]
  df$x <- relevel(factor(df$x), ref = "Normal")
  df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
  mod <- coxph(Surv(ftime, cenvar) ~x + GENDER + RACE_ETHNICITY + UC_SITE + ADI, data = df, weights = w) 
  coefs <- summary(mod)$coef
  out <- as.data.frame(coefs[grep("x", rownames(coefs)), c(2, 3, 6), drop = FALSE])
  out$Variable <- gsub("x", paste0(gsub("_cat", "", var, fixed = TRUE), ": "), rownames(out))
  out$Variable <- paste0(out$Variable, " vs. Normal")
  out <- dplyr::select(out, Variable, everything())
  return(out)
})
lab.results2 <- as.data.frame(do.call("rbind", lab.results))
lab.results2$adj.P.Val <- p.adjust(lab.results2$P, "bonferroni")
lab.results2$Sig <- ifelse(lab.results2$adj.P.Val < 0.05, "Significant", "")
display(lab.results2)

df_tab <- createDataFrame(lab.results2)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/labs_cause_specific.csv")
%md
## Variability
Variability models were adjusted for gender, race, ADI, and UC site, and most recent value of that lab.  Variability was estimated using variablity independent of the mean (VIM)

var <- "AST"
vim.results <- lapply(labs, function(var){
  df <- AD2
  df$vim <- df[[paste0(var, "_VIM_cat")]]
  df$recent <- df[[paste0(var, "_LATEST")]]
  # Model missingness
  df$nonmissing <- !is.na(df$vim)
  weights <- weightit(nonmissing ~ RACE_ETHNICITY + GENDER + AGE_LAST_OBS + ADI + UC_SITE, data = df, estimand = "ATE", method = "ps")
  df$w <- weights$weights
  df <- df[df$nonmissing,]
  df$vim <- relevel(factor(df$vim), ref = "q1")
  df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
  mod <- coxph(Surv(ftime, cenvar) ~vim + recent + GENDER + RACE_ETHNICITY + UC_SITE + ADI, data = df, weights = w) 
  coefs <- summary(mod)$coef
  out <- as.data.frame(coefs[grep("vim", rownames(coefs)), c(2, 3, 6), drop = FALSE])
  out$Variable <- gsub("vimq", paste0(var, ": Quintile "), rownames(out))
  out$Variable <- paste0(out$Variable, " vs. Quintile 1")
  out <- dplyr::select(out, Variable, everything())
  return(out)
})
vim.results2 <- as.data.frame(do.call("rbind", vim.results))
vim.results2$adj.P.Val <- p.adjust(vim.results2$P, "bonferroni")
vim.results2$Sig <- ifelse(vim.results2$adj.P.Val < 0.05, "Significant", "")
display(vim.results2)
df_tab <- createDataFrame(vim.results2)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/variability_cause_specific.csv")
%md
## Effects of BMI variability on AD by race
Had to exclude AI/AK native, Native Hawaiian/PI, and multirace subjects due to small numbers
table(AD2$BMI_VIM_cat, AD2$RACE_ETHNICITY)
df <- subset(AD2, !(RACE_ETHNICITY %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander", "Multirace")))
var <- "BMI_VIM_cat"
df$x <- df[[var]]
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~x*(RACE_ETHNICITY + GENDER) + UC_SITE + ADI, data = df) 
pairs0 <- pairs(emmeans(mod, ~x|RACE_ETHNICITY), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, RACE_ETHNICITY, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("RACE_ETHNICITY", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- subset(tmp, grepl("q1", contrast))
tmp <- dplyr::select(tmp, RACE_ETHNICITY, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_VIM_by_race")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("White", contrast1) & grepl("q1", contrast))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_VIM_race_interaction")
%md
## Effects of BMI variability on AD by gender
pairs0 <- pairs(emmeans(mod, ~x|GENDER), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, GENDER, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("GENDER", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- subset(tmp, grepl("q1", contrast))
tmp <- dplyr::select(tmp, GENDER, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_VIM_by_gender")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("q1", contrast))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_VIM_gender_interaction")
%md
## Effects of SBP variability on AD by race
Had to exclude AI/AK native and Native Hawaiian/PI subjects due to small numbers
df <- subset(AD2, !(RACE_ETHNICITY %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander")))
var <- "Sys_BP_3004249_VIM_cat"
df$x <- df[[var]]
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~x*(RACE_ETHNICITY + GENDER) + UC_SITE + ADI, data = df) 
pairs0 <- pairs(emmeans(mod, ~x|RACE_ETHNICITY), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, RACE_ETHNICITY, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("RACE_ETHNICITY", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- subset(tmp, grepl("q1", contrast))
tmp <- dplyr::select(tmp, RACE_ETHNICITY, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/SBP_VIM_by_race")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("White", contrast1) & grepl("q1", contrast))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/SBP_VIM_race_interaction")
%md
## Effects of SBP variability on AD by gender
pairs0 <- pairs(emmeans(mod, ~x|GENDER), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, GENDER, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("GENDER", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- subset(tmp, grepl("q1", contrast))
tmp <- dplyr::select(tmp, GENDER, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/SBP_VIM_by_gender")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("q1", contrast))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/SBP_VIM_gender_interaction")
%md
## Demographics
(GENDER, RACE_ETHNICITY, ADI, BMI_cat_CDC, BMI_cat_RACE)

Models for BMI were adjusted for gender, race, ADI, and UC site.  Other models were unadjusted.
df <- AD2
df$GENDER <- relevel(factor(df$GENDER), ref = "MALE")
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~GENDER, data = df) 
coefs <- summary(mod)$coef
out.GENDER <- as.data.frame(coefs[grep("GENDER", rownames(coefs)), c(2, 3, 5), drop = FALSE])
out.GENDER$Variable <- "Gender: Female vs. Male"
out.GENDER <- dplyr::select(out.GENDER, Variable, everything())
df <- AD2
df$RACE_ETHNICITY <- relevel(factor(df$RACE_ETHNICITY), ref = "White")
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~RACE_ETHNICITY, data = df) 
coefs <- summary(mod)$coef
out.race <- as.data.frame(coefs[grep("RACE_ETHNICITY", rownames(coefs)), c(2, 3, 5), drop = FALSE])
out.race$Variable <- gsub("RACE_ETHNICITY", "Race: ", rownames(out.race))
out.race$Variable <- paste0(out.race$Variable, " vs. White")
out.race <- dplyr::select(out.race, Variable, everything())
df <- AD2
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~ADI, data = df) 
coefs <- summary(mod)$coef
out.adi <- as.data.frame(coefs[grep("ADI", rownames(coefs)), c(2, 3, 5), drop = FALSE])
out.adi$Variable <- "ADI"
out.adi <- dplyr::select(out.adi, Variable, everything())
df <- subset(AD2, !is.na(AD2$BMI_MIN_cat_CDC))
df$BMI_MIN_cat_CDC <- relevel(factor(df$BMI_MIN_cat_CDC), ref = "Normal")
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~BMI_MIN_cat_CDC + GENDER + RACE_ETHNICITY + UC_SITE + ADI, data = df) 
coefs <- summary(mod)$coef
out.bmicdc.min <- as.data.frame(coefs[grep("BMI_MIN_cat_CDC", rownames(coefs)), c(2, 3, 5), drop = FALSE])
out.bmicdc.min$Variable <- gsub("BMI_MIN_cat_CDC", "BMI_MIN_cat_CDC: ", rownames(out.bmicdc.min))
out.bmicdc.min$Variable <- paste0(out.bmicdc.min$Variable, " vs. Normal")
out.bmicdc.min <- dplyr::select(out.bmicdc.min, Variable, everything())
df$BMI_MAX_cat_CDC <- relevel(factor(df$BMI_MAX_cat_CDC), ref = "Normal")
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~BMI_MAX_cat_CDC + GENDER + RACE_ETHNICITY + UC_SITE + ADI, data = df) 
coefs <- summary(mod)$coef
out.bmicdc.max <- as.data.frame(coefs[grep("BMI_MAX_cat_CDC", rownames(coefs)), c(2, 3, 5), drop = FALSE])
out.bmicdc.max$Variable <- gsub("BMI_MAX_cat_CDC", "BMI_MAX_cat_CDC: ", rownames(out.bmicdc.max))
out.bmicdc.max$Variable <- paste0(out.bmicdc.max$Variable, " vs. Normal")
out.bmicdc.max <- dplyr::select(out.bmicdc.max, Variable, everything())
out.bmicdc <- rbind(out.bmicdc.min, out.bmicdc.max)
df <- subset(AD2, !is.na(AD2$BMI_MIN_cat_RACE))
df$BMI_MIN_cat_RACE <- relevel(factor(df$BMI_MIN_cat_RACE), ref = "Normal")
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~BMI_MIN_cat_RACE + GENDER + RACE_ETHNICITY + UC_SITE + ADI, data = df) 
coefs <- summary(mod)$coef
out.bmirace.min <- as.data.frame(coefs[grep("BMI_MIN_cat_RACE", rownames(coefs)), c(2, 3, 5), drop = FALSE])
out.bmirace.min$Variable <- gsub("BMI_MIN_cat_RACE", "BMI_MIN_cat_RACE: ", rownames(out.bmirace.min))
out.bmirace.min$Variable <- paste0(out.bmirace.min$Variable, " vs. Normal")
out.bmirace.min <- dplyr::select(out.bmirace.min, Variable, everything())
df$BMI_MAX_cat_RACE <- relevel(factor(df$BMI_MAX_cat_RACE), ref = "Normal")
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~BMI_MAX_cat_RACE + GENDER + RACE_ETHNICITY + UC_SITE + ADI, data = df) 
coefs <- summary(mod)$coef
out.bmirace.max <- as.data.frame(coefs[grep("BMI_MAX_cat_RACE", rownames(coefs)), c(2, 3, 5), drop = FALSE])
out.bmirace.max$Variable <- gsub("BMI_MAX_cat_RACE", "BMI_MAX_cat_RACE: ", rownames(out.bmirace.max))
out.bmirace.max$Variable <- paste0(out.bmirace.max$Variable, " vs. Normal")
out.bmirace.max <- dplyr::select(out.bmirace.max, Variable, everything())
out.bmirace <- rbind(out.bmirace.min, out.bmirace.max)
out <- rbind(out.GENDER, out.race, out.adi, out.bmicdc, out.bmirace)
display(out)
df_tab <- createDataFrame(out)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/demog_cause_specific.csv")
%md
## Effect of lowest BMI (CDC) on AD by race
Had to exclude AI/AK native and Native Hawaiian/PI subjects due to low numbers
df <- subset(AD2, !(RACE_ETHNICITY %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander")))
var <- "BMI_MIN_cat_CDC"
df$x <- df[[var]]
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~x*(RACE_ETHNICITY + GENDER) + UC_SITE + ADI, data = df) 
pairs0 <- pairs(emmeans(mod, ~x|RACE_ETHNICITY), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, RACE_ETHNICITY, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("RACE_ETHNICITY", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, RACE_ETHNICITY, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MIN_cat_CDC_by_race")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("White", contrast1))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MIN_cat_CDC_race_interaction")
%md
## Effect of lowest BMI (CDC) on AD by gender
pairs0 <- pairs(emmeans(mod, ~x|GENDER), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, GENDER, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("GENDER", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, GENDER, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MIN_cat_CDC_by_gender")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MIN_cat_CDC_gender_interaction")
%md
## Effect of highest BMI (CDC) on AD by race
Subjects who were underweight at their highest BMI were excluded due to low numbers, as were AI/AK native and native Hawaiian/PI subject
df <- subset(AD2, !(RACE_ETHNICITY %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander")) & BMI_MAX_cat_CDC != "Underweight")
var <- "BMI_MAX_cat_CDC"
df$x <- factor(df[[var]], levels = c("Normal", "Overweight", "Obese"))
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~x*(RACE_ETHNICITY + GENDER) + UC_SITE + ADI, data = df) 
pairs0 <- pairs(emmeans(mod, ~x|RACE_ETHNICITY), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, RACE_ETHNICITY, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("RACE_ETHNICITY", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, RACE_ETHNICITY, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MAX_cat_CDC_by_race")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("White", contrast1))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MAX_cat_CDC_race_interaction")
%md
## Effects of highest BMI (CDC) on AD by gender
pairs0 <- pairs(emmeans(mod, ~x|GENDER), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, GENDER, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("GENDER", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, GENDER, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MAX_cat_CDC_by_gender")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MAX_cat_CDC_gender_interaction")
%md
## Effects of lowest BMI (Race-based) on AD by race
df <- subset(AD2, !(RACE_ETHNICITY %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander")))
var <- "BMI_MIN_cat_RACE"
df$x <- df[[var]]
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~x*(RACE_ETHNICITY + GENDER) + UC_SITE + ADI, data = df) 
pairs0 <- pairs(emmeans(mod, ~x|RACE_ETHNICITY), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, RACE_ETHNICITY, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("RACE_ETHNICITY", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, RACE_ETHNICITY, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MIN_cat_RACE_by_race")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("White", contrast1))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MIN_cat_RACE_race_interaction")
%md 
## Effect of lowest BMI (race-based) on AD by gender
pairs0 <- pairs(emmeans(mod, ~x|GENDER), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, GENDER, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("GENDER", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, GENDER, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MIN_cat_RACE_by_gender")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MIN_cat_RACE_gender_interaction")
%md
## Effect of highest BMI (race-based) on AD by race
df <- subset(AD2, !(RACE_ETHNICITY %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander")) & BMI_MAX_cat_RACE != "Underweight")
var <- "BMI_MAX_cat_RACE"
df$x <- df[[var]]
df$x <- factor(df$x, levels = c("Normal", "Overweight", "Obese"))
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~x*(RACE_ETHNICITY + GENDER) + UC_SITE + ADI, data = df) 
pairs0 <- pairs(emmeans(mod, ~x|RACE_ETHNICITY), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, RACE_ETHNICITY, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("RACE_ETHNICITY", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, RACE_ETHNICITY, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MAX_cat_RACE_by_race")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("White", contrast1))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MAX_cat_RACE_race_interaction")
%md
## Effect of maximum BMI (race-based) on AD by gender
pairs0 <- pairs(emmeans(mod, ~x|GENDER), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, GENDER, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("GENDER", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, GENDER, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MAX_cat_RACE_by_gender")
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/BMI_MAX_cat_RACE_gender_interaction")
%md
## Table of frequency of all diagnoses make these categorical
tabdat <- AD2
vars.use <- diagnoses
tmp <- lapply(diagnoses, function(var){
  tabdat[[var]] <<- ifelse(!is.na(tabdat[[var]]), "Yes", "No")
})
tabdat$Status <- factor(tabdat$fstatus, levels = c(1, 2, 0))
levels(tabdat$Status) <- c("AD", "Died Without AD", "Lost to Followup")
tabdat$HOT_FLASHES <- ifelse(tabdat$GENDER == "Male", NA, tabdat$HOT_FLASHES)
tabdat$HOT_FLASHES_MP <- ifelse(tabdat$GENDER == "Male", NA, tabdat$HOT_FLASHES_MP)
form <- paste0("~", paste(vars.use, collapse = "+"), "|Status")
tab <- as.data.frame(table1(as.formula(form), data = tabdat))
display(tab)
df_tab <- createDataFrame(tab)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/diagnoses_by_group")
%md
## Cross-tabulation of ALCOHOL by NONINFECTIOUS_HEPATITIS
tabdat <- AD2
tabdat$ALCOHOL <- ifelse(!is.na(tabdat$ALCOHOL), "Yes", "No")
tabdat$NONINFECTIOUS_HEPATITIS <- ifelse(!is.na(tabdat$NONINFECTIOUS_HEPATITIS), "Yes", "No")
tab <- as.data.frame(table(Alcohol = tabdat$ALCOHOL, `Noninfectious Hepatitis` = tabdat$NONINFECTIOUS_HEPATITIS))
display(tab)
df_tab <- createDataFrame(tab)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/alcohol_by_noninfectious_hepatitis")
%md
## Gender by race
df <- subset(AD2, !(RACE_ETHNICITY %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander")))
df$GENDER <- relevel(factor(df$GENDER), ref = "MALE")
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~ GENDER*RACE_ETHNICITY, data = df) 
pairs0 <- pairs(emmeans(mod, ~GENDER|RACE_ETHNICITY), reverse = TRUE, adjust = "none")
pairs.df <- as.data.frame(pairs0)
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, RACE_ETHNICITY, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("RACE_ETHNICITY", "contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, RACE_ETHNICITY, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, Comparison = contrast, `Hazard Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)

df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/gender_by_race")
%md
## Differences between races in effect of gender
pairs <- pairs(pairs0, by = "contrast", adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, contrast1, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast", "contrast1"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast1, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("White", contrast1))
tmp <- dplyr::rename(tmp, `Comparison 1` = contrast1, `Comparison 2` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/gender_race_interaction")
%md
## ADI by race and gender
df <- subset(AD2, !(RACE_ETHNICITY %in% c("American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander")))
df$GENDER <- relevel(factor(df$GENDER), ref = "MALE")
df$cenvar <- ifelse(df$fstatus == 1, 1, 0)
mod <- coxph(Surv(ftime, cenvar) ~ ADI*(GENDER + RACE_ETHNICITY), data = df) 
pairs0 <- emtrends(mod, "RACE_ETHNICITY", var = "ADI")
pairs0
pairs.df <- as.data.frame(test(pairs0, adjust = "none"))
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, RACE_ETHNICITY, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("RACE_ETHNICITY"))
tmp$estimate <- round(exp(tmp$ADI.trend), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, RACE_ETHNICITY, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, `Hazard Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/ADI_by_race")
%md
## Differences between races in effect of ADI
pairs <- pairs(pairs0, adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- subset(tmp, grepl("White", contrast))
tmp <- dplyr::rename(tmp, `Comparison` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/ADI_race_interaction")
%md
## ADI by gender
pairs0 <- emtrends(mod, "GENDER", var = "ADI")
pairs0
pairs.df <- as.data.frame(test(pairs0, adjust = "none"))
ci.df <- as.data.frame(confint(pairs0))
ci.df <- dplyr::select(ci.df, GENDER, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("GENDER"))
tmp$estimate <- round(exp(tmp$ADI.trend), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, GENDER, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, `Hazard Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/ADI_by_gender")
%md 
## Differences between genders in effect of ADI
pairs <- pairs(pairs0, adjust = "none")
pairs.df <- as.data.frame(pairs)
ci.df <- as.data.frame(confint(pairs))
ci.df <- dplyr::select(ci.df, contrast, asymp.LCL, asymp.UCL)
tmp <- dplyr::left_join(pairs.df, ci.df, by = c("contrast"))
tmp$estimate <- round(exp(tmp$estimate), 5)
tmp$p.value <- format.pval(tmp$p.value, digits = 6, eps = 0.00001)
tmp$conf.low <- round(exp(tmp$asymp.LCL), 5)
tmp$conf.high <- round(exp(tmp$asymp.UCL), 5)
tmp <- dplyr::select(tmp, contrast, estimate, conf.low, conf.high, p.value) 
tmp <- dplyr::rename(tmp, `Comparison` = contrast, `Odds Ratio` = estimate, `P-Value` = p.value, `95% CI Lower Bound` = conf.low, `95% CI Upper Bound` = conf.high)
display(tmp)
df_tab <- createDataFrame(tmp)
createOrReplaceTempView(df_tab, "df_tab")
%python
df_tab = spark.sql("select * from df_tab")
df_tab.show()
df_tab.coalesce(1).write.format("csv").option("header","true").mode("OverWrite").csv("dbfs:/mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/ADI_gender_interaction")
%fs
ls /mnt/ucd-bpdurbin-workspace/Yvonne_Wan/cause_specific_2023-11-13/