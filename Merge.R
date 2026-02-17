#Creation of the dataset with data from OpenMeteo and RegioneLombardia

#Libraries----
library(readr)
library(lubridate)
library(dplyr)

#Import Data----
#Import data of daily pollution concentrations (PM2.5 and PM10) from the stations in Milan
datavalore <- read_csv("Dati_sensori_aria_2010-2017_20251109.csv")
#Import data about the stations
stazioni <-read_csv("stazioni_milano.csv")
names(stazioni)[names(stazioni) == "IdSensore"] <- "idSensore"

#Pollution Data----
datavalore$Data <- as.POSIXct(datavalore$Data, format = "%d/%m/%Y %H:%M:%S")
datavalore$Data <- as.Date(datavalore$Data)
datavalore$week <- isoweek(datavalore$Data)
datavalore$year <- year(datavalore$Data)

#We aggregate data by week
#for each week we take the average of the daily pollution if we have at least 30% of valid values
#otherwise we set the weekly mean to NA
weekly_mean <- datavalore %>%
  group_by(idSensore, year, week) %>%
  summarise(
    avg_val = if (sum(Valore>=0)/7 < 0.3) NA_real_ else mean(Valore[Valore>=0])
  )

#Weather Data----
#Import daily weather data from OpenMeteo
open_meteo_all <- read_csv("open-meteo-all.csv")

open_meteo_all$time <- as.Date(open_meteo_all$time)
open_meteo_all$week <- isoweek(open_meteo_all$time)
open_meteo_all$year <- year(open_meteo_all$time)

open_meteo_all$week[open_meteo_all$week ==53 & open_meteo_all$year==2010] = 0

#aggregate by week
#weekly min(max) temperature = min(max) of daily min(max) temperature
#weekly rain = sum of daily rain
#weekly wind speed = max of daily wind speed

weekly_summary_cov <- open_meteo_all %>%
  group_by(location_id, year, week) %>%
  summarise(
    temp_max_max = max(`temperature_2m_max (°C)`, na.rm = TRUE),

    temp_min_min = min(`temperature_2m_min (°C)`, na.rm = TRUE),

    rain_sum_total = sum(`rain_sum (mm)`, na.rm = TRUE),
    
    wind_max = max(`wind_speed_10m_max (km/h)`, na.rm = TRUE),

  )

# we change the location_id to match id_sensore
weekly_summary_cov$location_id[weekly_summary_cov$location_id ==0] =10320
weekly_summary_cov$location_id[weekly_summary_cov$location_id ==1] =20429
weekly_summary_cov$location_id[weekly_summary_cov$location_id ==2] =17122
weekly_summary_cov$location_id[weekly_summary_cov$location_id ==3] =10273
weekly_summary_cov$location_id[weekly_summary_cov$location_id ==4] =6956
weekly_summary_cov$location_id[weekly_summary_cov$location_id ==5] =10283
weekly_summary_cov$location_id[weekly_summary_cov$location_id ==6] =9874
weekly_summary_cov$location_id[weekly_summary_cov$location_id ==7] =6905
weekly_summary_cov$location_id[weekly_summary_cov$location_id ==8] =20529


names(weekly_summary_cov)[names(weekly_summary_cov) == "location_id"] <- "idSensore"

#Merge pollution and weather Data----

merge_weekly_cov2 = merge(weekly_summary_cov, weekly_mean, c("idSensore", "week", "year"))

merge_weekly_cov2$pollutant[merge_weekly_cov2$idSensore == 10320] = "PM10"
merge_weekly_cov2$pollutant[merge_weekly_cov2$idSensore == 20429] = "PM10"
merge_weekly_cov2$pollutant[merge_weekly_cov2$idSensore == 17122] = "PM2.5"
merge_weekly_cov2$pollutant[merge_weekly_cov2$idSensore == 10273] = "PM10"
merge_weekly_cov2$pollutant[merge_weekly_cov2$idSensore == 6956] = "PM10"
merge_weekly_cov2$pollutant[merge_weekly_cov2$idSensore == 10283] = "PM2.5"
merge_weekly_cov2$pollutant[merge_weekly_cov2$idSensore == 9874] = "PM2.5"
merge_weekly_cov2$pollutant[merge_weekly_cov2$idSensore == 6905] = "PM10"
merge_weekly_cov2$pollutant[merge_weekly_cov2$idSensore == 20529] = "PM2.5"

#Save Dataset----
write_csv(merge_weekly_cov2, "merged_data.csv")



