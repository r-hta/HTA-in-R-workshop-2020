# Script to rearrange bugs_loghr and bugs_baseline to have same order and use same
# column names as event_names

event_names <- c("Stroke", "MI", "Bleed", "ICH", "Death", "SE", "TIA")
treatment_names <- c("Coumarin", "Apixaban", "Dabigatran", "Edoxaban")


# Load the baseline/coumarin log hazard and treatment log hazard ratios estimated by a 
# network meta-analysis in OpenBUGS
bugs_log_hr <- read.csv(file = "bugs_loghr.csv")
bugs_log_hazard <- read.csv(file = "bugs_baseline.csv")

bugs_log_hr_temp <- matrix(NA, nrow = dim(bugs_log_hr)[1], ncol = dim(bugs_log_hr)[2])
bugs_log_hazard_temp <- matrix(NA, nrow = dim(bugs_log_hazard)[1], ncol = dim(bugs_log_hazard)[2])

# Log hazard event_name
colnames(bugs_log_hazard_temp) <- paste("Log hazard", event_names)
# Log HR Treatment event_name
colnames(bugs_log_hr_temp) <- paste("Log HR", rep(treatment_names[-1], each = length(event_names)),
                                        rep(event_names, 3))


for(event_name in event_names) {
  bugs_log_hazard_temp[, paste("Log hazard", event_name)] <- 
    bugs_log_hazard[, grepl(event_name, colnames(bugs_log_hazard), ignore.case = TRUE)]
  for(treatment_name in treatment_names[-1]) {
    bugs_log_hr_temp[, paste("Log HR", treatment_name, event_name)] <-
      bugs_log_hr[, grepl(event_name, colnames(bugs_log_hr), ignore.case = TRUE) &
                    grepl(treatment_name, colnames(bugs_log_hr), ignore.case = TRUE)]
    
  }
}

# Now overwrite the previous data
write.csv(bugs_log_hr_temp, file = "bugs_loghr.csv")
write.csv(bugs_log_hazard_temp, file = "bugs_baseline.csv")

# Needed to edit in Excel to remove empty first column.
