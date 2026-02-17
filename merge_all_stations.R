#Creation of the dataset for spatial stan model

#Libraries----
library(lubridate)
library(dplyr)
library(ggplot2)

#Upload pollution data----
all_stations <- read.csv("Stazioni_PM10_Lombardia.csv")
names(all_stations)[names(all_stations) == "IdSensore"] <- "idSensore"
datavalore <- read.csv("data_all_stations.csv")


datavalore$Data <- as.POSIXct(datavalore$Data, format = "%d/%m/%Y %H:%M:%S")
datavalore$Data <- as.Date(datavalore$Data)
datavalore$week <- isoweek(datavalore$Data)
datavalore$year <- year(datavalore$Data)


#weekly average of pollution data
weekly_mean <- datavalore %>%
  group_by(idSensore, year, week) %>%
  summarise(
    avg_val = if (sum(Valore>=0)/7 < 0.3) NA_real_ else mean(Valore[Valore>=0])
  )

weekly_mean$week <- as.numeric(weekly_mean$week)
weekly_mean <- weekly_mean %>%
  mutate(
    data = as.Date(paste(year, week, 1, sep = "-"), format = "%Y-%W-%u")
  )

#plot time series of each sensor to see which ones are missing a lot of values
i <- 1
#go through the list of idSensore
current_sensor <- stazioni$idSensore[i]
#plot
ggplot(weekly_mean[weekly_mean$idSensore == current_sensor, ], aes(x = data, y = avg_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(title = "weekly pollutant", x = "week", y = "pollutant") +
  theme_minimal()
#increase index
i <- i+1
#if the current sensor is missing a lot of data take it out of the dataset
weekly_mean <- weekly_mean[which(weekly_mean$idSensore != current_sensor), ]

unique(weekly_mean$idSensore)
#we are left with 51 sensors.

#sensors: 6908  6910  6912  6915  6917  6918  6919  6921  6925  6926  6927  6931  6935  6947
#  6951  6952  6953  6954  6955  6956  6958  6959  6962  9739  9831  9835  9872  9890
#  9927  9961  9963  9965  9976 10027 10043 10072 10083 10103 10173 10273 10290 10320
# 10339 10343 10352 10483 10484 10501 10548 20197 20565


write_csv(weekly_mean, "weekly_all_stations.csv")

weekly_all_stations <- read.csv("weekly_all_stations.csv")
sum(is.na(weekly_all_stations$avg_val))

sensors_id <- c(6908,  6910,  6912,  6915,  6917,  6918,  6919,  6921,  6925,  6926,  6927,  6931,  6935,  6947,
                6951,  6952,  6953,  6954,  6955,  6956,  6958,  6959,  6962,  9739,  9831,  9835,  9872,  9890,
                9927,  9961,  9963,  9965,  9976, 10027, 10043, 10072, 10083, 10103, 10173, 10273, 10290, 10320,
                10339, 10343, 10352, 10483, 10484, 10501, 10548, 20197, 20565)

stations <- all_stations[which(all_stations$IdSensore %in% sensors_id), ]
write_csv(stations, "stations.csv")

#Stations coordinates for Open-Meteo----
stations <- read.csv("stations.csv")
cat(
  paste(stations$lat, stations$lng, sep = ", "),
  sep = "\n"
)
#45.8169745, 8.82024911
#45.62735705, 9.02640138
#45.478347, 9.231667
#45.64971, 9.60123469
#45.8085738, 9.2217792
#45.01688067, 11.07609728
#45.54006456, 10.22281989
#45.08199189, 9.70079278
#45.470501, 9.19746075
#45.73083728, 9.12573936
#45.46242074, 8.88021388
#45.15175186, 10.78141717
#46.135486, 9.56624
#45.58289429, 8.84216546
#45.04650286, 11.18094889
#45.64962671, 10.20509022
#45.36631075, 9.70395139
#45.660992, 9.160225
#46.16797103, 9.87014581
#45.48363242, 9.32736196
#45.61924753, 8.75697656
#45.69043647, 9.48426111
#45.51306214, 10.19190606
#45.85019631, 9.39558125
#45.14254358, 10.04384767
#45.69568669, 9.66126061
#44.99955481, 9.00844961
#45.69104328, 9.64366092
#45.51933533, 9.59201917
#45.14526864, 10.80333489
#45.81504217, 9.06697886
#45.30278678, 9.49528589
#45.16057658, 10.79557331
#45.63387247, 9.55609572
#45.46334886, 9.19532511
#45.23349689, 9.66625758
#45.574117, 9.263808
#45.51571417, 10.33624003
#45.64605175, 10.38929672
#45.18633544, 9.14667703
#45.15937678, 9.69815314
#45.28194289, 8.75433614
#46.46952397, 10.37543275
#45.87460256, 10.17736553
#45.84221272, 9.35166247
#45.49822794, 9.55624181
#45.52000944, 9.51223153
#45.41277633, 10.68336214
#45.10277464, 8.90418742
#45.52655611, 8.73650211
#45.86375372, 9.39995481

#Upload weather data----
weather_locations <- read.csv("open-meteo-locations.csv")
stations <- cbind(stations, weather_locations$location_id)

#merge weekly data of the pollutants with the location of each sensor
weekly_all_stations <- merge(
  weekly_all_stations,
  stations[, c("idSensore", "lat", "lng", "Utm_Nord", "UTM_Est", "weather_locations$location_id")],
  by = "idSensore",
  all.x = TRUE
)
colnames(weekly_all_stations) <- c("idSensore","year","week","avg_val",                      
                                   "data","lat","lng","Utm_Nord","UTM_Est","location_id")
write_csv(weekly_all_stations, "weekly_all_stations.csv")

#weekly computation for weather data
open_meteo_all <- read_csv("open-meteo-region.csv")


open_meteo_all$time <- as.Date(open_meteo_all$time)
open_meteo_all$week <- isoweek(open_meteo_all$time)
open_meteo_all$year <- year(open_meteo_all$time)

open_meteo_all$week[open_meteo_all$week ==53 & open_meteo_all$year==2010] =0

weekly_summary_cov <- open_meteo_all %>%
  group_by(location_id, year, week) %>%
  summarise(
    temp_max_max = max(`temperature_2m_max (°C)`, na.rm = TRUE),
    
    temp_min_min = min(`temperature_2m_min (°C)`, na.rm = TRUE),
    
    rain_sum_total = sum(`rain_sum (mm)`, na.rm = TRUE),
    
    wind_max = max(`wind_speed_10m_max (km/h)`, na.rm = TRUE),
    
  )

data_region <- merge(weekly_all_stations, weekly_summary_cov, c("location_id", "week", "year"))

#Saving the Data---- 
write_csv(data_region, "data_region.csv")
