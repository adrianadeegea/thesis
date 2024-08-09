rm()

# Installing packages
install.packages("readr")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("lubridate")
install.packages("tidyr")
install.packages("stringr")
install.packages("broom")
install.packages("gridExtra")
install.packages("hms")
install.packages("nlme")
install.packages("tidyverse")
install.packages("car")
install.packages("ggpubr")
install.packages("reshape2")
install.packages("ggdist")
install.packages("Matrix")
install.packages("MuMIn")
install.packages("multcomp")
install.packages("knitr")

# Loading libraries
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(stringr)
library(broom)
library(gridExtra)
library(hms)
library(nlme)
library(tidyverse)
library(car)
library(ggpubr)
library(reshape2)
library(ggdist)
library(Matrix)
library(MuMIn)
library(multcomp)
library(knitr)

# First Part: Daily Analysis of Soil Respiration Flux 

# Importing dataframe
data1 <- May_Daily_Plot_C

# Keeping only useful rows and columns 
data1 <- data1[data1$`Tag(M3)` == 'M5', !(names(data1) %in% c("Tsen", "H2O", "O2", "Error", "Aux_V", "PAR", "Pressure", "Flow", "Rec_No", "Date", "Tag(M3)"))]

# Renaming Plot_No column
data1 <- data1 %>% rename(Plot = Plot_No)

# Converting Plot numbers to categorical variables
data1$Plot <- as.factor(data1$Plot)

# Converting the Time column to hms format
data1 <- data1 %>%
  mutate(Time = as_hms(Time))

# Extracting the hour from the Time column
data1 <- data1 %>%
  mutate(time = paste0(ifelse(substr(Time, 1, 2) == "0" | substr(Time, 1, 2) == "00", substr(Time, 1, 1), substr(Time, 1, 2)), ":00:00"),
         time = case_when(
           substr(Time, 1, 2) == "08" ~ "09:00:00",
           substr(Time, 1, 2) == "18" ~ "19:00:00",
           substr(Time, 1, 2) == "20" ~ "21:00:00",
           TRUE ~ time ))

# Converting CO2 from ppm to mg/m3
data1$mg_per_m3_per_s <- data1$CO2 * 44.01 / 24.45

# Converting from mg/m3 to microgram/m3
data1$ug_per_m3_per_s <- data1$mg_per_m3_per_s * 1000

# Converting from microgram/m3 to microgram/cm3
data1$ug_per_cm3_per_s <- data1$ug_per_m3_per_s / 1e6

# Creating a function to calculate the slope
calculate_slope <- function(df) {
  lm(ug_per_cm3_per_s ~ Time, data = df) %>%
    tidy() %>%
    filter(term == "Time") %>%
    pull(estimate)}

# Grouping by time and Plot and calculating the CO2 slope and Tsoil Tair and Msoil averages
result <- data1 %>%
  group_by(time, Plot) %>%
  summarize(
    slope = calculate_slope(cur_data()),
    Tsoil = mean(Tsoil, na.rm = TRUE),
    Tair = mean(Tair, na.rm = TRUE),
    Msoil = mean(Msoil, na.rm = TRUE)
  ) %>%
  ungroup()

# Creating a new dataframe with the results
data2 <- result

# Filtering the volumes dataframe to include only rows where month is "May"
may_volumes <- volumes %>%
  filter(Month == "May") %>%
  select(Plot, volume, Month)
may_volumes$Plot <- as.factor(may_volumes$Plot)

# Merging the volume column to data2
data2 <- data2 %>%
  left_join(may_volumes %>% select(Plot, volume), by = "Plot")

# Setting Area 
data2 <- data2 %>%
  mutate(area = 78)

# Calculating flux in umol/cm2/s
data2 <- data2 %>%
  mutate(flux_ug_cm3 = slope * volume / area)

# Converting slopes from microgram/cm³/s to micromol/m³/s
data2$flux_umol_m2_sec <- data2$flux_ug_cm3 / (44.01/10000)

# Assigning treatment type
data2 <- data2 %>%
  mutate(Treatment = case_when(
    Plot == 1 ~ "N",
    Plot == 2 ~ "NIB",
    Plot == 21 ~ "C",
    Plot == 22 ~ "S",
    Plot == 23 ~ "SI",
    Plot == 24 ~ "NI",
    Plot == 25 ~ "SB",
    Plot == 26 ~ "NB",
    Plot == 27 ~ "SIB"))

# Creating a new dataframe for ANOVA
data3 <- data2

# Converting Time column to character
data3$time <- as.character(data3$time)
data3 <- data3 %>% rename(Time = time)
data3 <- data3 %>% rename(Flux = flux_umol_m2_sec)

# Creating columns based on Treatment
data3$basalt_application <- NA
data3$tree_type <- NA
data3$inoculation <- NA
data3$basalt_application <- ifelse(data3$Treatment %in% c("NB", "SB", "NIB", "SIB"), "Yes", "No")
data3$basalt_application[data3$Treatment == "C"] <- "pretreatment"
data3$tree_type <- "None"
data3$tree_type[data3$Treatment %in% c("N", "NB", "NIB", "NI")] <- "Native Broadleaf"
data3$tree_type[data3$Treatment %in% c("SB", "S", "SIB", "SI")] <- "Spruce"
data3$tree_type[data3$Treatment == "C"] <- "pretreatment"
data3$inoculation <- ifelse(data3$Treatment %in% c("NIB", "SIB", "SI", "NI"), "Yes", "No")
data3$inoculation[data3$Treatment == "C"] <- "pretreatment"

# Converting columns Time, basalt_application, and tree_type to factors
data3$basalt_application <- as.factor(data3$basalt_application)
data3$tree_type <- as.factor(data3$tree_type)
data3$Time <- as.factor(data3$Time)

# Converting columns Flux, Tsoil, and Msoil to numericals
data3$Flux <- as.numeric(data3$Flux)
data3$Tsoil <- as.numeric(data3$Tsoil)
data3$Msoil <- as.numeric(data3$Msoil)

# Checking if the data meets the ANOVA assumptions
OP <- par(mfrow=c(1,2))
plot(lm(Flux ~ basalt_application * tree_type * inoculation + Tsoil + Msoil + Time, data = data3), 1:2)
par(OP)

# Converting Flux to logaritmic function to meet ANOVA assumptions
data3$log_Flux <- log(data3$Flux)

# Checking if the new data meets the ANOVA assumptions
OP <- par(mfrow=c(1,2))
plot(lm(log_Flux ~ basalt_application * tree_type * inoculation + Tsoil + Msoil + Time, data = data3), 1:2)
par(OP)

# Fitting the ANOVA model
daily_model <- aov(log_Flux ~ basalt_application * tree_type * inoculation + Tsoil + Msoil + Error(Time), data = data3)

# Summarizing the results
summary(daily_model)

# Post-Hoc for Time (TukeyHSD)
TukeyHSD(aov(log_Flux ~ Time, data = data3))

# Calculating daily averages
average_Flux <- data3 %>%
  group_by(Treatment) %>%
  summarize(average_Flux = mean(Flux, na.rm = TRUE)) %>%
  arrange(desc(average_Flux))

# Printing daily averages
print(average_Flux)

# Creating a plot for Treatment C
dataC <- data3 %>%
  filter(Treatment %in% c("C"))

dataC_plot <- ggplot(dataC, aes(x = Time, y = Flux, color = Treatment, group = Treatment)) +
  geom_line(size = 0.5) +
  geom_point(size = 2) +
  labs(title = "Control", x = "Time", y = "Rs Flux (µmol/m²/s)") +
  scale_color_manual(values = c("C" = "orange"), 
                     labels = c("C" = "C")) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "#D3D3D3"),
        panel.grid.minor = element_line(color = "#D3D3D3"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.title.x = element_text(vjust = -0.5), 
        axis.title.y = element_text(vjust = 2, margin = margin(r = 20))) + 
  scale_y_continuous(limits = c(1.5, 5.5), breaks = seq(1.5, 5.5, by = 0.5))

# Creating a plot for Treatment N
dataN <- data3 %>%
  filter(Treatment %in% c("N", "NB", "NI", "NIB"))

dataN_plot <- ggplot(dataN, aes(x = Time, y = Flux, color = Treatment, group = Treatment)) +
  geom_line(size = 0.5) +
  geom_point(size = 2) +
  labs(title = "Native Broadleaf", x = "Time", y = "Rs Flux (µmol/m²/s)") +
  scale_color_manual(values = c("N" = "lightblue", "NB" = "#1f78b4", "NI" = "#6a3d9a", "NIB" = "darkblue"), 
                     labels = c("N"="N", "NB"="NB", "NI"="NI", "NIB"="NIB")) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "#D3D3D3"),
        panel.grid.minor = element_line(color = "#D3D3D3"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.title.x = element_text(vjust = -0.5), 
        axis.title.y = element_text(vjust = 2, margin = margin(r = 20))) + 
  scale_y_continuous(limits = c(1.5, 5.5), breaks = seq(1.5, 5.5, by = 0.5))

# Creating a plot for Treatment S
dataS <- data3 %>%
  filter(Treatment %in% c("S", "SB", "SI", "SIB"))

dataS_plot <- ggplot(dataS, aes(x = Time, y = Flux, color = Treatment, group = Treatment)) +
  geom_line(size = 0.5) +
  geom_point(size = 2) +
  labs(title = "Spruce", x = "Time", y = "Rs Flux (µmol/m²/s)") +
  scale_color_manual(values = c("S" = "lightgreen", "SB" = "#66c2a5", "SI" = "#2ca25f", "SIB" = "darkgreen"), 
                     labels = c("S"="S", "SB"="SB", "SI"="SI", "SIB"="SIB")) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "#D3D3D3"),
        panel.grid.minor = element_line(color = "#D3D3D3"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.title.x = element_text(vjust = -0.5), 
        axis.title.y = element_text(vjust = 2, margin = margin(r = 20))) + 
  scale_y_continuous(limits = c(1.5, 5.5), breaks = seq(1.5, 5.5, by = 0.5))

# Creating the main title
main_title <- textGrob("Daily Variation of Soil Respiration (Rs) Flux ", x = 0.5, hjust = 0.5, gp = gpar(fontsize = 17, fontface = "bold"))

# Combining the plots with the main title
combined_plot <- grid.arrange(main_title, 
                              dataC_plot, 
                              dataN_plot, 
                              dataS_plot, 
                              nrow = 4, 
                              heights = c(0.1, 1, 1, 1)) 

# Printing the combined plot
print(combined_plot)


----------------------------------------------------------------------------------------------------------
  
  
# Second Part: Seasonal Analysis of Soil Respiration Flux 

# Importing dataframes
august <- X1_August_2023_Clean_1 
october <- X2_October_2023_clean_1
january <- X3_January_2024_clean_1 
may <- X4_May_2024_clean_1_csv 

# Cleaning data
august <- august %>% filter(`Tag(M3)` == "M5")
october <- october %>% filter(`Tag(M3)` == "M5")
may <- may %>% filter(`Tag(M3)` == "M5")
january <- january %>% filter(`Tag(M3)` == "M5")

columns_to_delete <- c("Tsen", "H2O", "O2", "Error", "Aux_V", "PAR", "Pressure", "Flow", "Rec_No", "Date", "Tag(M3)")
august <- august %>% select(-all_of(columns_to_delete))
october <- october %>% select(-all_of(columns_to_delete))
may <- may %>% select(-all_of(columns_to_delete))
january <- january %>% select(-all_of(columns_to_delete))

# Merging dataframes into Month1
may <- may %>% mutate(Month = "May")
august <- august %>% mutate(Month = "August")
october <- october %>% mutate(Month = "October")
january <- january %>% mutate(Month = "January")
months1 <- bind_rows(may, august, october, january)

# Renaming Plot_No column
months1 <- months1 %>% rename(Plot = Plot_No)

# Converting Plot numbers to categorical variables
months1$Plot <- as.factor(months1$Plot)

# Converting the Time column to hms format
months1 <- months1 %>%
  mutate(Time = as_hms(Time))

# Converting CO2 from ppm to mg/m3
months1$mg_per_m3_per_s <- months1$CO2 * 44.01 / 24.45

# Converting from mg/m3 to microgram/m3
months1$ug_per_m3_per_s <- months1$mg_per_m3_per_s * 1000

# Converting from microgram/m3 to microgram/cm3
months1$ug_per_cm3_per_s <- months1$ug_per_m3_per_s / 1e6

# Calculating slopes per Plot within each Month
slopes <- months1 %>%
  group_by(Month, Plot) %>%
  do({
    model <- lm(ug_per_cm3_per_s ~ Time, data = .)
    slope <- coef(model)[2]  
    data.frame(Slope = slope)
  }) %>%
  ungroup()

# Calculating averages per Plot within each Month
averages <- months1 %>%
  group_by(Month, Plot) %>%
  summarise(
    Tsoil = mean(Tsoil, na.rm = TRUE),
    Tair = mean(Tair, na.rm = TRUE),
    Msoil = mean(Msoil, na.rm = TRUE)
  ) %>%
  ungroup()

# Merging the averages and slopes dataframes to create months2
months2 <- slopes %>%
  left_join(averages, by = c("Month", "Plot"))

# Adding volumes in cm³
volumes_selected <- volumes %>%
  select(Plot, Month, volume)
volumes_selected$Plot <- as.factor(volumes_selected$Plot)

# Merging volumes_selected with months2 based on 'Plot' and 'Month'
months2 <- left_join(months2, volumes_selected, by = c("Plot", "Month"))

# Setting Area in cm²
months2 <- months2 %>%
  mutate(area = 78)

# Calculating flux
months2 <- months2 %>%
  mutate(Flux = Slope * volume / area)

# Converting slopes from microgram/cm³/s to micromol/m³/s
months2$flux_umol_m2 <- months2$Flux / (44.01/10000)

# Creating Block column  
months2 <- months2 %>%
  mutate(Block = case_when(
    Plot %in% c(1, 2, 21, 22, 23, 24, 25, 26, 27) ~ "C",
    Plot %in% c(3, 4, 5, 6, 7, 8, 9, 10, 11) ~ "A",
    Plot %in% c(46, 47, 48, 49, 50, 51, 52, 53, 54) ~ "F",
    Plot %in% c(64, 65, 66, 67, 68, 69, 70, 71, 72) ~ "H",
    TRUE ~ NA_character_  ))

# Assigning a treatment to each Plot
months2 <- months2 %>%
  mutate(Treatment = case_when(
    Plot == 1 ~ "N",
    Plot == 2 ~ "NIB",
    Plot == 3 ~ "N",
    Plot == 4 ~ "NIB",
    Plot == 5 ~ "SIB",
    Plot == 6 ~ "SI",
    Plot == 7 ~ "NB",
    Plot == 8 ~ "SB",
    Plot == 9 ~ "C",
    Plot == 10 ~ "S",
    Plot == 11 ~ "NI",
    Plot == 21 ~ "C",
    Plot == 22 ~ "S",
    Plot == 23 ~ "SI",
    Plot == 24 ~ "NI",
    Plot == 25  ~ "SB",
    Plot == 26 ~ "NB",
    Plot == 27 ~ "SIB",
    Plot == 46 ~ "NB",
    Plot == 47 ~ "N",
    Plot == 48 ~ "C",
    Plot == 49 ~ "NIB",
    Plot == 50 ~ "SIB",
    Plot == 51 ~ "S",
    Plot == 52 ~ "NI",
    Plot == 53 ~ "SB",
    Plot == 54 ~ "SI",
    Plot == 64 ~ "S",
    Plot == 65 ~ "NB",
    Plot == 66 ~ "N",
    Plot == 67 ~ "C",
    Plot == 68 ~ "SIB",
    Plot == 69 ~ "NIB",
    Plot == 70 ~ "SI",
    Plot == 71 ~ "NI",
    Plot == 72 ~ "SB"))

# Creating a new dataframe for ANOVA
months3 <- months2

# Creating columns based on treatments
months3$basalt_application <- NA
months3$tree_type <- NA
months3$inoculation <- NA
months3$basalt_application <- ifelse(months3$Treatment %in% c("NB", "SB", "NIB", "SIB"), "Yes", "No")
months3$basalt_application[months3$Treatment == "C"] <- "pretreatment"
months3$tree_type <- "None"
months3$tree_type[months3$Treatment %in% c("N", "NB", "NIB", "NI")] <- "Native Broadleaf"
months3$tree_type[months3$Treatment %in% c("SB", "S", "SIB", "SI")] <- "Spruce"
months3$tree_type[months3$Treatment == "C"] <- "pretreatment"
months3$inoculation <- ifelse(months3$Treatment %in% c("NIB", "SIB", "SI", "NI"), "Yes", "No")
months3$inoculation[months3$Treatment == "C"] <- "pretreatment"

# Converting columns Month, basalt_application, and tree_type to factors
months3$basalt_application <- as.factor(months3$basalt_application)
months3$tree_type <- as.factor(months3$tree_type)
months3$Month <- as.factor(months3$Month)

# Converting columns log_Flux, Tsoil, and Msoil to numericals
months3$Flux <- as.numeric(months3$Flux)
months3$Tsoil <- as.numeric(months3$Tsoil)
months3$Msoil <- as.numeric(months3$Msoil)

# Checking if the results meet the ANOVA assumptions
OP <- par(mfrow=c(1,2))
plot(lm(Flux ~ basalt_application * tree_type * inoculation + Tsoil + Msoil+ Month, data = months3), 1:2)
par(OP)

# Converting Flux to logaritmic function to meet ANOVA assumptions
months3$log_Flux <- log(months3$Flux + 1) 

# Removing outliers
months3 <- months3[-78, , drop = FALSE]
months3 <- months3[-18, , drop = FALSE]
months3 <- months3[-2, , drop = FALSE]
months3 <- months3[-1, , drop = FALSE]

# Checking if the results meet the ANOVA assumptions
OP <- par(mfrow=c(1,2))
plot(lm(log_Flux ~ basalt_application * tree_type * inoculation + Tsoil + Msoil+ Month, data = months3), 1:2)
par(OP)

# Fitting the ANOVA model
monthly_model <- aov(log_Flux ~ basalt_application * tree_type * inoculation + Tsoil + Msoil+ Error(Block/Month), data = months3)

# Summarizing the results
summary(monthly_model)

# Calculating annual averages and standard errors
average_Flux2 <- months3 %>%
  group_by(Treatment) %>%
  summarize(
    average_Flux = mean(flux_umol_m2, na.rm = TRUE),
    sd_Flux = sd(flux_umol_m2, na.rm = TRUE),
    n = n(),
    se_Flux = sd_Flux / sqrt(n)
  ) %>%
  arrange(desc(average_Flux))

# Printing annual averages with standard errors
print(average_Flux2)

# Calculating R-squared for the entire model
SST_Block <- 4.297e-05 + 9.692e-05  
SST_Block_Month <- 0.001175 + 0.003514 + 0.000015 + 0.001456  
SST_Within <- 0.0001059 + 0.0000122 + 0.0000381 + 0.0000005 + 0.0000018 + 0.0000026 + 0.0000001 + 0.0000018 + 0.0008987  

# Calculating Total Sum of Squares
SST <- SST_Block + SST_Block_Month + SST_Within

# Calculating Residual Sum of Squares
SSres_Block_Month <- 0.001456  
SSres_Within <- 0.0008987 
SSres <- SSres_Block_Month + SSres_Within

# Calculating Regression Sum of Squares
SSreg <- SST - SSres

# Calculating R-squared for the entire model
R_squared_overall <- SSreg / SST

# Calculating R-squared only for significant variables
SSM <- 0.003514 + 0.0001059 + 0.0000381  

# Calculating R-squared for significant variables
R_squared_significant <- SSM / SST

# Calculating Percentage of Variability Explained by Each Significant Variable
percentage_Tsoil_Block_Month <- 0.003514 / SST * 100
percentage_basalt_application_Within <- 0.0001059 / SST * 100
percentage_inoculation_Within <- 0.0000381 / SST * 100

# Printing all Results
cat("R-squared for the entire model:", R_squared_overall, "\n")
cat("R-squared for significant variables:", R_squared_significant, "\n")
cat("Percentage of Variability Explained by Tsoil in Block:Month:", percentage_Tsoil_Block_Month, "%\n")
cat("Percentage of Variability Explained by basalt_application in Within:", percentage_basalt_application_Within, "%\n")
cat("Percentage of Variability Explained by inoculation in Within:", percentage_inoculation_Within, "%\n")

# Specifying Month order for the Seasonal Plots
months3$Month <- factor(months3$Month, levels = c("August", "October", "January", "May"))

# Creating plot for treatments with Native Broadleaf
filtered_dataN <- months3 %>%
  filter(Treatment %in% c("C", "N", "NB", "NI", "NIB"))

summary_dataN <- filtered_dataN %>%
  group_by(Month, Treatment) %>%
  summarise(
    mean_flux = mean(flux_umol_m2, na.rm = TRUE),
    se_flux = sd(flux_umol_m2, na.rm = TRUE) / sqrt(n()))

plotsN <- summary_dataN %>%
  ggplot(aes(x = Treatment, y = mean_flux, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.6) +
  geom_errorbar(aes(ymin = mean_flux - se_flux, ymax = mean_flux + se_flux), 
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_grid(. ~ Month, scales = "free_y", switch = "y") + 
  labs(x = "Treatment",
       y = "Rs Flux (\u00B5mol m\u207B\u00B2)",
       fill = "Treatment") + 
  scale_fill_manual(values = c("C" = "orange", 
                               "N" = "lightblue", 
                               "NB" = "#377EB8", 
                               "NI" = "blue", 
                               "NIB" = "darkblue")) +
  theme_minimal() +
  theme(legend.position = "right", 
        panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12, face = "bold"), 
        plot.title = element_blank(), 
        axis.title.x = element_text(vjust = -1.5), 
        axis.title.y = element_text(vjust = 2, margin = margin(r = 20)), 
        axis.text.x = element_text(hjust = 0.5), 
        axis.text.y = element_text(size = 10), 
        panel.grid.major = element_line(color = "grey70"), 
        panel.grid.minor = element_line(color = "grey85")) 

# Creating plot for treatments with Spruce
filtered_dataS <- months3 %>%
  filter(Treatment %in% c("C", "S", "SB", "SI", "SIB"))

summary_dataS <- filtered_dataS %>%
  group_by(Month, Treatment) %>%
  summarise(
    mean_flux = mean(flux_umol_m2, na.rm = TRUE),
    se_flux = sd(flux_umol_m2, na.rm = TRUE) / sqrt(n()))

plotsS <- summary_dataS %>%
  ggplot(aes(x = Treatment, y = mean_flux, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.6) +
  geom_errorbar(aes(ymin = mean_flux - se_flux, ymax = mean_flux + se_flux), 
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_grid(. ~ Month, scales = "free_y", switch = "y") + 
  labs(x = "Treatment",
       y = "Rs Flux (\u00B5mol m\u207B\u00B2)",
       fill = "Treatment") + 
  scale_fill_manual(values = c("C" = "orange", 
                               "S" = "lightgreen", 
                               "SB" = "#4DAF4A", 
                               "SI" = "green", 
                               "SIB" = "darkgreen")) +
  theme_minimal() +
  theme(legend.position = "right", 
        panel.spacing = unit(2, "lines"), 
        strip.text.x = element_text(size = 12, face = "bold"), 
        plot.title = element_blank(), 
        axis.title.x = element_text(vjust = -1.5), 
        axis.title.y = element_text(vjust = 2, margin = margin(r = 20)), 
        axis.text.x = element_text(hjust = 0.5), 
        axis.text.y = element_text(size = 10), 
        panel.grid.major = element_line(color = "grey70"), 
        panel.grid.minor = element_line(color = "grey85")) 

# Creating row titles
row_title_N <- textGrob("Native Broadleaf", x = 0.05, hjust = 0, gp = gpar(fontsize = 14))
row_title_S <- textGrob("Spruce", x = 0.05, hjust = 0, gp = gpar(fontsize = 14))

# Creating the main title
main_title <- textGrob("Seasonal Variation of Soil Respiration (Rs) Flux", x = 0.5, hjust = 0.5, gp = gpar(fontsize = 17, fontface = "bold"))

# Combining the plots with the main title and row titles
combined_plot <- grid.arrange(main_title, 
                              arrangeGrob(row_title_N, plotsN, ncol = 1, heights = c(0.1, 1)),
                              arrangeGrob(row_title_S, plotsS, ncol = 1, heights = c(0.1, 1)),
                              ncol = 1, 
                              heights = c(0.1, 1.1, 1.1)) 

# Printing the combined plot
print(combined_plot)
