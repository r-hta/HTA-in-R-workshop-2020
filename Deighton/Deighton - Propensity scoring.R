#SET OPTIONS AND LOAD PACKAGES-----
options(scipen=999)
options(digits = 3)
set.seed(123)

library("install.load")
packages <- c("dplyr", #data handling
              "weights") #weighted t-test
install_load(packages)
rm(packages)

#CREATE DATASETS-----
#Create dataset for Drug A
bl_age <- rnorm(n = 500, mean = 36, sd = 15)
bl_bmi <- rnorm(n = 500, mean = 24, sd = 3)
bl_creatinine <- rnorm(n = 500, mean = 330, sd = 230)
bl_platelets <- rnorm(n = 500, mean = 115, sd = 60)
bl_ldh <- rnorm(n = 500, mean = 480, sd = 520)
bl_egfr <- rnorm(n = 500, mean = 25, sd = 15)

ep_creatinine <- rnorm(n = 500, mean = 150, sd = 75)
ep_platelets <- rnorm(n = 500, mean = 235, sd = 65)
ep_ldh <- rnorm(n = 500, mean = 180, sd = 35)
ep_egfr <- rnorm(n = 500, mean = 55, sd = 30)

df_a <- data.frame(cbind(bl_age, bl_bmi, bl_creatinine, bl_platelets, bl_ldh, bl_egfr,
                   ep_creatinine, ep_platelets, ep_ldh, ep_egfr))

#Create dataset for Drug B
bl_age <- rnorm(n = 500, mean = 37, sd = 16)
bl_bmi <- rnorm(n = 500, mean = 23, sd = 4)
bl_creatinine <- rnorm(n = 500, mean = 400, sd = 240)
bl_platelets <- rnorm(n = 500, mean = 120, sd = 85)
bl_ldh <- rnorm(n = 500, mean = 710, sd = 580)
bl_egfr <- rnorm(n = 500, mean = 24, sd = 22)

ep_creatinine <- rnorm(n = 500, mean = 158, sd = 100)
ep_platelets <- rnorm(n = 500, mean = 245, sd = 75)
ep_ldh <- rnorm(n = 500, mean = 195, sd = 60)
ep_egfr <- rnorm(n = 500, mean = 70, sd = 35)

df_b <- data.frame(cbind(bl_age, bl_bmi, bl_creatinine, bl_platelets, bl_ldh, bl_egfr,
                     ep_creatinine, ep_platelets, ep_ldh, ep_egfr))

#Combine datasets, remove old data and create cuts of the dataset
df <- rbind(df_a, df_b)
df$drug <- factor(c(rep("A", 500), rep("B", 500)))
rm(list = setdiff(ls(), "df"))

for (i in seq(from = 0, to = 400, by = 100)) {
  assign(
    paste0("df", (i/100)+1), 
    df[c((1 + i):(100 + i), (nrow(df) - i):(nrow(df) - (i + 99))), ]
  )
}

#START ANALYSIS-----
for (a in 1:5) {
  
  #Select dataframe
  if (a == 1) {analysis_data <- df1}
  if (a == 2) {analysis_data <- df2}
  if (a == 3) {analysis_data <- df3}
  if (a == 4) {analysis_data <- df4}
  if (a == 5) {analysis_data <- df5}
  
  #Perform propensity scoring-----
  prop_model <- glm(drug ~  bl_creatinine + bl_platelets + bl_ldh + bl_egfr,
                    data = analysis_data,
                    family = binomial(link="logit"),
                    control = list(maxit = 50))
  analysis_data$pr_score <- predict(prop_model, type="response")
  p <- nrow(analysis_data[analysis_data$drug == "B",]) / nrow(analysis_data)
  analysis_data$SW <- ifelse(analysis_data$drug == "B", p / analysis_data$pr_score, (1 - p) / (1 - analysis_data$pr_score))
  
  #Tabulate baseline data-----
  vars_bl <- c("bl_age", "bl_bmi", "bl_creatinine", "bl_platelets", "bl_ldh", "bl_egfr")
  for (b in (1 : length(vars_bl))) {
    name_b <- vars_bl[b]
    x <- analysis_data[analysis_data$drug == "A", name_b]
    y <- analysis_data[analysis_data$drug == "B", name_b]
    test <- t.test(y, x)
    x_desc <- data.frame("Mean" = mean(x, na.rm=T))
    y_desc <- data.frame("Mean" = mean(y, na.rm=T))
    xy_desc <- rbind(x_desc,y_desc)
    xy_desc <- xy_desc %>% t() %>% 
      as.data.frame() %>% 
      `colnames<-`(c("Drug A", "Drug B")) %>% 
      mutate("p-value" = test$p.value,
             "variable" = name_b) %>% 
      select("variable", "Drug A", "Drug B", "p-value")
    
    table_name <- paste("BL_df", a, name_b, sep = "_")
    assign(table_name, xy_desc)
  }
  
  #Tabulate baseline data after applying stabilised weights-----
  vars_bl <- c("bl_age", "bl_bmi", "bl_creatinine", "bl_platelets", "bl_ldh", "bl_egfr")
  for (b in (1 : length(vars_bl))) {
    name_b <- vars_bl[b]
    x <- analysis_data[analysis_data$drug == "A", c(name_b, "SW")]
    y <- analysis_data[analysis_data$drug == "B", c(name_b, "SW")]
    test <- wtd.t.test(x = y[, name_b], 
                       y = x[, name_b], 
                       weight = y[, "SW"], 
                       weighty = x[, "SW"],
                       mean1 = T, 
                       samedata = F)
    x_desc <- data.frame("Mean" = wtd.mean(x = x[, name_b], weights = x[, "SW"], na.rm = T))
    y_desc <- data.frame("Mean" = wtd.mean(x = y[, name_b], weights = y[, "SW"], na.rm = T))
    xy_desc <- rbind(x_desc,y_desc)
    xy_desc <- xy_desc %>% t() %>% 
      as.data.frame() %>% 
      `colnames<-`(c("Drug A", "Drug B")) %>% 
      mutate("p-value" = test$coefficients[3],
             "variable" = name_b) %>% 
      select("variable", "Drug A", "Drug B", "p-value")
    
    table_name <- paste("BL_SW_df", a, name_b, sep = "_")
    assign(table_name, xy_desc)
  }

  list_bl_unadjusted <- mget(ls(pattern="BL_df")) #put dataframes containing "BL_df" in a list
  list_bl_SW <- mget(ls(pattern="BL_SW_df")) #put dataframes containing "BL_SW_df" in a list
  BL_unadjusted <- do.call(what = rbind, args = list_bl_unadjusted) #combine dfs and exclude list descriptor from list of dfs to bind  
  BL_SW <- do.call(what = rbind, args = list_bl_SW) #combine dfs and exclude list descriptor from list of dfs to bind  

  df_name <- paste("BL_unadj", a, sep="_")
  assign(df_name, BL_unadjusted)
  df_name <- paste("BL_SW", a, sep="_")
  assign(df_name, BL_SW)

  rm(list=(ls(pattern="BL_df"))) #remove datasets to prevent overlap between loops
  rm(list=(ls(pattern="BL_SW_df"))) #remove datasets to prevent overlap between loops
  rm(BL_SW, BL_unadjusted) #remove datasets to prevent overlap between loops
 
  #Tabulate endpoint data-----
  vars_ep <- c("ep_creatinine", "ep_platelets", "ep_ldh", "ep_egfr")
  for (b in (1 : length(vars_ep))) {
    name_b <- vars_ep[b]
    x <- analysis_data[analysis_data$drug == "A", name_b]
    y <- analysis_data[analysis_data$drug == "B", name_b]
    test <- t.test(y, x)
    x_desc <- data.frame("Mean" = mean(x, na.rm=T))
    y_desc <- data.frame("Mean" = mean(y, na.rm=T))
    xy_desc <- rbind(x_desc,y_desc)
    xy_desc <- xy_desc %>% t() %>% 
      as.data.frame() %>% 
      `colnames<-`(c("Drug A", "Drug B")) %>% 
      mutate("p-value" = test$p.value,
             "variable" = name_b) %>% 
      select("variable", "Drug A", "Drug B", "p-value")
    
    table_name <- paste("EP_df", a, name_b, sep = "_")
    assign(table_name, xy_desc)
  }
  
  #Tabulate baseline data after applying stabilised weights-----
  vars_ep <- c("ep_creatinine", "ep_platelets", "ep_ldh", "ep_egfr")
  for (b in (1 : length(vars_ep))) {
    name_b <- vars_ep[b]
    x <- analysis_data[analysis_data$drug == "A", c(name_b, "SW")]
    y <- analysis_data[analysis_data$drug == "B", c(name_b, "SW")]
    test <- wtd.t.test(x = y[, name_b], 
                       y = x[, name_b], 
                       weight = y[, "SW"], 
                       weighty = x[, "SW"],
                       mean1 = T, 
                       samedata = F)
    x_desc <- data.frame("Mean" = wtd.mean(x = x[, name_b], weights = x[, "SW"], na.rm = T)) 
    y_desc <- data.frame("Mean" = wtd.mean(x = y[, name_b], weights = y[, "SW"], na.rm = T)) 
    xy_desc <- rbind(x_desc,y_desc)
    xy_desc <- xy_desc %>% t() %>% 
      as.data.frame() %>% 
      `colnames<-`(c("Drug A", "Drug B")) %>% 
      mutate("p-value" = as.numeric(test$coefficients[3]),
             "variable" = name_b) %>% 
      select("variable", "Drug A", "Drug B", "p-value")
    
    table_name <- paste("EP_SW_df", a, name_b, sep = "_")
    assign(table_name, xy_desc)
  }
  
  list_ep_unadjusted <- mget(ls(pattern="EP_df")) #put dataframes containing "BL_df" in a list
  list_ep_SW <- mget(ls(pattern="EP_SW_df")) #put dataframes containing "BL_SW_df" in a list
  EP_unadjusted <- do.call(what = rbind, args = list_ep_unadjusted) #combine dfs and exclude list descriptor from list of dfs to bind  
  EP_SW <- do.call(what = rbind, args = list_ep_SW) #combine dfs and exclude list descriptor from list of dfs to bind  
  
  df_name <- paste("EP_unadj", a, sep="_")
  assign(df_name, EP_unadjusted)
  df_name <- paste("EP_SW", a, sep="_")
  assign(df_name, EP_SW)
  
  rm(list=(ls(pattern="EP_df"))) #remove datasets to prevent overlap between loops
  rm(list=(ls(pattern="EP_SW_df"))) #remove datasets to prevent overlap between loops
  rm(EP_SW, EP_unadjusted) #remove datasets to prevent overlap between loops
} 

#Assess datasets-----
rm(list=ls()[! ls() %in% c("BL_SW_1","BL_SW_2","BL_SW_3","BL_SW_4","BL_SW_5",
                           "BL_unadj_1","BL_unadj_2","BL_unadj_3","BL_unadj_4","BL_unadj_5",
                           "EP_SW_1","EP_SW_2","EP_SW_3","EP_SW_4","EP_SW_5",
                           "EP_unadj_1","EP_unadj_2","EP_unadj_3","EP_unadj_4","EP_unadj_5"
                           )])

print(BL_unadj_4, row.names = F)
print(EP_unadj_4, row.names = F)
print(BL_SW_4, row.names = F)
print(EP_SW_4, row.names = F)
