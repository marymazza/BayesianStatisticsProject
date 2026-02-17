#Exploratory analysis

#Libraries----
library(lubridate)
library(dplyr)
library(ggplot2)
library(corrplot)

#Upload Weekly merged Data----
merged_data <- read.csv("merged_data.csv")
merged_data$week <- as.numeric(merged_data$week)

merged_data <- merged_data %>%
  mutate(
    data = as.Date(paste(year, week, 1, sep = "-"), format = "%Y-%W-%u")
  )


#Plotting the time series----

#Temperature range
x11()
ggplot(merged_data, aes(x = data)) +
  geom_ribbon(aes(ymin = temp_min_min, ymax = temp_max_max), fill = "grey", alpha = 0.3) +
  geom_line(aes(y = temp_max_max), color = "red", linewidth = 1) +
  geom_line(aes(y = temp_min_min), color = "blue", linewidth = 1) +
  labs(
    title = "weekly temperature range",
    x = "week",
    y = "Temperature (°C)"
  ) +
  theme_minimal()

#Rain
ggplot(merged_data, aes(x = data, y = rain_sum_total)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  labs(title = "weekly rain", x = "week", y = "rain") +
  theme_minimal()

#Wind speed
ggplot(merged_data, aes(x = data, y = wind_max)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly wind", x = "week", y = "wind") +
  theme_minimal()

#Pollution
ggplot(merged_data[merged_data$idSensore == 10273, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant sensor 10273", x = "week", y = "pollutant") +
  theme_minimal()

ggplot(merged_data[merged_data$idSensore == 10320, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant sensor 10320", x = "week", y = "pollutant") +
  theme_minimal()

ggplot(merged_data[merged_data$idSensore == 17122, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant sensor 17122", x = "week", y = "pollutant") +
  theme_minimal()
#NOTE: sensor 17122 only recorded data from the end of 2012 to the end of 2017

ggplot(merged_data[merged_data$idSensore == 6956, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant sensor 6956", x = "week", y = "pollutant") +
  theme_minimal()

ggplot(merged_data[merged_data$idSensore == 10283, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant sensor 10283", x = "week", y = "pollutant") +
  theme_minimal()

ggplot(merged_data[merged_data$idSensore == 9874, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant sensor 9874", x = "week", y = "pollutant") +
  theme_minimal()
#NOTE: sensor 9874 has all NA values

ggplot(merged_data[merged_data$idSensore == 6905, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant sensor 6905", x = "week", y = "pollutant") +
  theme_minimal()
#NOTE: sensor 6905 has all NA values

ggplot(merged_data[merged_data$idSensore == 20529, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant sensor 20529", x = "week", y = "pollutant") +
  theme_minimal()
ggplot(merged_data[merged_data$idSensore == 20429, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant sensor 20429", x = "week", y = "pollutant") +
  theme_minimal()
#NOTE: sensor 20429 and 20529 didn't record any data

#pairs plot
pairs(merged_data[,-c(1,2,3,9,10)], col = as.factor(merged_data$idSensore))

#Deleting useless data----

#deleting all sensors that don't have all data and adding a column for inside/outside area c indicator

new_merged_data <- merged_data[which(merged_data$idSensore %in% c("10273", "6956", "10283")), ]
new_merged_data$inside <- ifelse(new_merged_data$idSensore == "6956", 1, 0)

#since we have only one sensor (10283) that measures PM2.5 we can't use it for Causal Impact
#we work on PM10: consider sensor 6956 for data inside area c and sensor 10273 for data outside area c

rm(data)
data <- cbind (new_merged_data[which(new_merged_data$idSensore == "6956"), 10],
               new_merged_data[which(new_merged_data$idSensore == "6956"), 8],
               new_merged_data[which(new_merged_data$idSensore == "10273"), 4:8])
colnames(data) <- c("date","avg_val_inside", "T_max", "T_min", "rain", "wind", "avg_val_outside")
data <- data[!is.na(data$date), ]

#Correlation plot----
cor_matrix <- cor(data[,-c(1,8)], use = "pairwise.complete.obs")
x11()
corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  tl.col = "black",
  addCoef.col = "black",
  number.cex = 1.8,
  col = colorRampPalette(c("#B2182B", "white", "#2166AC"))(200)
)

#high correlation between levels of PM10 inside and outside area c
#as expected wind and rain decrease levels of PM10

#Save Dataset----
write_csv(data, "causalImpact_data.csv")


