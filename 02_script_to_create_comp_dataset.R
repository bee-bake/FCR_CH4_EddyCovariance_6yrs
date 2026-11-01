#Script to create a comprehensive dataset of daily methane and CO2 fluxes with other environmental and lake drivers
#created by: Bibek Kandel on 2025 Feb 02
#Reviewed by ABP on 2025 July 23
#Updated by BK on 2025 July 25

#Load libraries
pacman::p_load(tidyverse,ggplot2,lubridate,rLakeAnalyzer, vroom)

###########
# This chunk of code will read in high-frequency meteorological data file from EDI and break them into 30 minute intervals to match them to half-hourly Eddy covariance data. 
#Pull EDI file
met_EDI <- vroom(
  "https://pasta.lternet.edu/package/data/eml/edi/389/10/d3f3d2fa40c41fdcd505ae49b2fdcf8b?key=yltMpS4UEIk12AvB9L7OL5uRiG0",
  col_select = c("DateTime","PAR_umolm2s_Average","Rain_Total_mm","WindSpeed_Average_m_s",
    "AirTemp_C_Average"),
  col_types = list(
    DateTime = col_datetime()
  )
) |>
  filter(DateTime >= ymd_hms("2020-04-04 00:15:00")) |>
  mutate(DateTime = with_tz(DateTime, "America/New_York"))

#Current met data
#Pull from LTREB site
#https://www.ltreb-reservoirs.org/data/
met_L1 <- vroom(
  "https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data-qaqc/FCRmet_L1.csv",
  col_select = c("DateTime","PAR_umolm2s_Average","Rain_Total_mm","WindSpeed_Average_m_s",
                 "AirTemp_C_Average"),
  col_types = list(
    DateTime = col_datetime()
  )
) |>
  filter(DateTime <= ymd_hms("2026-04-04 00:15:00")) |>
  mutate(DateTime = with_tz(DateTime, "America/New_York"))

#merge
met_all <- bind_rows(met_EDI, met_L1)

# Select data every 30 minutes from Jan 2020 to end of met data
met_all$Breaks <- cut(met_all$DateTime,breaks = "30 mins",right=FALSE)
met_all$Breaks <- ymd_hms(as.character(met_all$Breaks))

# Average met data to the 30 min mark (excluding Total Rain and Total PAR)
met_30 <- met_all |> 
  group_by(Breaks) |> 
  summarise_all(mean,na.rm=TRUE)

# average to the daily time step
met_daily <- met_30 |> 
  mutate(date=as.Date(DateTime)) |> 
  group_by(date) |> 
  summarise_all(mean,na.rm=TRUE) |> 
  select(-Breaks, -DateTime)
############


############
#This chunk of code will load in quality controlled but not gap filled Eddy covariance data from the github repository, filter days with more than 20 half hourly fluxes and sum them to calculate the daily fluxes.  
ec <- read_csv("C:/Users/bibek/Desktop/Bibek_files/5_yrs_FCR_flux_analysis/Eddy_fcr_footprint_full.csv", show_col_types = FALSE) |>
  select(-...1)

# ec <- read_csv("https://raw.github.com/bee-bake/5_yrs_FCR_flux_analysis/main/Eddy_fcr_footprint_full.csv", show_col_types = FALSE) |> 
#   select(-...1)   # Removes the unnecessary first column

ec2 <- ec |> 
  mutate(DateTime = as.POSIXct(datetime, "%Y-%m-%dT%H:%M:%SZ", tz = "EST"))

#In the following steps, we only included the days where there are more than or equal to 20 half hourly fluxes to calculate daily methane fluxes.
ec_selected <- ec2 |> 
  select(DateTime, ch4_flux_umolm2s, co2_flux_umolm2s) |> 
  dplyr::filter(DateTime >= "2020-04-04" & DateTime <= "2026-04-04") |> 
  mutate(date = as.Date(DateTime))

# Count the number of data points per day for each variable
ch4_co2_counts <- ec_selected |> 
  select(date, ch4_flux_umolm2s, co2_flux_umolm2s) |> 
  group_by(date) |> 
  summarise(
    ch4_count = sum(!is.na(ch4_flux_umolm2s)),
    co2_count = sum(!is.na(co2_flux_umolm2s)),
    .groups = "drop"
  )

#Inner join counts with ec_selected
ec_more_than_20 <- inner_join(ec_selected, ch4_co2_counts, by = "date")

#calculate the daily flux values
ec_daily_ch4 <- ec_more_than_20 |> 
  dplyr::filter(ch4_count >= 20) |> #this step filters days with half hourly flux values more than or equal to 15
  group_by(date) |> 
  summarise(daily_CH4 = mean(ch4_flux_umolm2s, na.rm = TRUE))

#calculate the daily flux values
ec_daily_co2 <- ec_more_than_20 |> 
  dplyr::filter(co2_count >= 20 ) |> #this step filters days with half hourly flux values more than or equal to 15
  group_by(date) |> 
  summarise(daily_CO2 = mean(co2_flux_umolm2s, na.rm = TRUE))

ec_daily <- full_join(ec_daily_co2, ec_daily_ch4, by = "date")

############


############
#This chunk of code will load in dissolved CH4 and CO2 data from EDI and calculate the mean daily values from two replicates at each depth.

##Filter the date and depth needed and convert to fluxes
ghg_EDI <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/551/10/ee8e65a2a380c9fb63bbf7f4a542c895?key=yltMpS4UEIk12AvB9L7OL5uRiG0") |> 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) |>  
  filter(Reservoir == "FCR" & Site == 50) |> 
  filter(DateTime >= "2020-04-04") |>  
  mutate(DateTime = round_date(DateTime, "30 mins")) 


ghg_L1 <- read.csv(
  "https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataNotYetUploadedToEDI/Raw_GHG/L1_manual_GHG.csv",
  stringsAsFactors = FALSE
) |> 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) |>  
  filter(Reservoir == "FCR" & Site == 50) |> 
  filter(DateTime <= "2026-04-04") |>  
  mutate(DateTime = round_date(DateTime, "30 mins"))

ghg_all <- rbind(ghg_EDI, ghg_L1)

## Plot to check
## Methane
ggplot(ghg_EDI,mapping=aes(x=DateTime,y=CH4_umolL))+
  geom_point()+
  geom_line()+
  theme_classic(base_size = 15)

## Carbon dioxide
ggplot(ghg_EDI,mapping=aes(x=DateTime,y=CO2_umolL))+
  geom_point()+
  geom_line()+
  theme_classic(base_size = 15)

## Pivot wider
ghg_fluxes_FCR <- ghg_EDI |> 
  dplyr::mutate(date = as.Date(DateTime)) |> 
  dplyr::select(date, Rep, CH4_umolL, CO2_umolL, Depth_m) |> 
  dplyr::group_by(date, Depth_m) |> 
  dplyr::summarize(mean_CH4 = mean(CH4_umolL, na.rm = TRUE),
                   mean_CO2 = mean(CO2_umolL, na.rm = TRUE),
                   .groups = "drop")

############


############
#This chunk of code will load in catwalk data from EDI, take out all midnight observations and evaluate half-hourly average values to match eddy covariance fluxes

catwalk_edi <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/271/10/814580ebec0385c66f0a0a97c38e9136?key=yltMpS4UEIk12AvB9L7OL5uRiG0")|> 
# takes out all midnight observations
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) |> 
  filter(DateTime >= "2020-04-04")

catwalk_L1 <- read_csv("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data-qaqc/fcre-waterquality_L1.csv") %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M:%S", tz="EST"))) |> 
  filter(DateTime <= "2026-04-04")

#Bind L1 file
catwalk_edi <- bind_rows(catwalk_edi, catwalk_L1)

# Calculate average half-hourly values
# Select data every 30 minutes from Jan 2020 to end of met data
catwalk_edi$Breaks <- as.character(cut(catwalk_edi$DateTime,breaks = "30 mins",right=FALSE))
catwalk_edi$Breaks <- lubridate::parse_date_time(catwalk_edi$Breaks, orders = c('ymd HMS','ymd')) 
# Average to half-hourly catwalk measurements
  catwalk_30 <- catwalk_edi |>
  select(DateTime,RDO_mgL_5_adjusted, RDO_mgL_9_adjusted, EXODO_mgL_1, Breaks, ThermistorTemp_C_1,
            ThermistorTemp_C_5, ThermistorTemp_C_9, EXOfDOM_QSU_1, EXOChla_ugL_1, EXOTemp_C_1) |> 
  group_by(Breaks) |>
  summarise_all(mean,na.rm=TRUE)
# Remove existing Datetime and use breaks as new half-hourly time in EST  
catwalk_30 <- catwalk_30 |> 
  select(Breaks,RDO_mgL_5_adjusted, RDO_mgL_9_adjusted, EXODO_mgL_1, ThermistorTemp_C_1,
             ThermistorTemp_C_5, ThermistorTemp_C_9, EXOfDOM_QSU_1, EXOChla_ugL_1, EXOTemp_C_1) |> 
  mutate(Breaks = force_tz(Breaks,"EST"))
#rename 'breaks' to 'DateTime'
  names(catwalk_30)[names(catwalk_30) == 'Breaks'] <- 'DateTime'
water_temp_9m <- catwalk_edi |>
  mutate(date = as.Date(DateTime),
             Month = lubridate::month(date)) |>
  select(date, EXOTemp_C_1, ThermistorTemp_C_5, ThermistorTemp_C_9, EXOfDOM_QSU_1, EXOChla_ugL_1, EXODOsat_percent_1, RDOsat_percent_5_adjusted, RDOsat_percent_9_adjusted,
             EXODO_mgL_1, RDO_mgL_5_adjusted, RDO_mgL_9_adjusted) |>
  na.omit() |>
  group_by(date) |>
  summarise(temp_1.6m = mean(EXOTemp_C_1),
                temp_5m = mean(ThermistorTemp_C_5),
                temp_9m = mean(ThermistorTemp_C_9),
                DOsat_surface = mean(EXODOsat_percent_1),
                DOsat_5m = mean(RDOsat_percent_5_adjusted),
                DOsat_9m = mean(RDOsat_percent_9_adjusted),
                DOmgL_surface = mean(EXODO_mgL_1),
                DOmgL_5m = mean(RDO_mgL_5_adjusted),
                DOmgL_9m = mean(RDO_mgL_9_adjusted),
                mean_Chla = mean(EXOChla_ugL_1),
                mean_fDOM = mean(EXOfDOM_QSU_1)) |>
  select(date, temp_1.6m, temp_5m, temp_9m, mean_fDOM, mean_Chla, DOsat_surface, DOsat_5m, DOsat_9m, DOmgL_surface, DOmgL_5m, DOmgL_9m)
#############


#############
#This chunk of code will prepare thermocline depth and buoyancy dataset using LakeAnalyzer package in R and catwalk temperature data.

# Format catwalk temp data for use in LakeAnalyzer. 
catwalk_temp <- catwalk_edi |> 
  select(DateTime,ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_3,
         ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, ThermistorTemp_C_7, ThermistorTemp_C_8,
         ThermistorTemp_C_9) |> 
  group_by(by30 = cut(DateTime, "30 min")) |> 
  summarise_all(mean,na.rm=TRUE) |>
  select(-DateTime) |>
  rename(datetime = by30, wtr_0.1 = ThermistorTemp_C_surface, wtr_1.0 = ThermistorTemp_C_1, 
         wtr_2.0 = ThermistorTemp_C_2, wtr_3.0 = ThermistorTemp_C_3, wtr_4.0 = ThermistorTemp_C_4,
         wtr_5.0 = ThermistorTemp_C_5, wtr_6.0 = ThermistorTemp_C_6, wtr_7.0 = ThermistorTemp_C_7,
         wtr_8.0 = ThermistorTemp_C_8, wtr_9.0 = ThermistorTemp_C_9)


# Using this function because it does dates with time and then just ymd for midnight
catwalk_temp$datetime <- lubridate::parse_date_time(catwalk_temp$datetime, orders = c('ymd HMS','ymd'))

fcr_temp = data.frame(catwalk_temp)


# Find the depth of the thermocline
thermo_depth =rLakeAnalyzer::ts.thermo.depth(fcr_temp,seasonal=TRUE)  

# Replace NaNs with NA
thermo_depth <- thermo_depth |>
  mutate(thermo.depth = ifelse(thermo.depth == "NaN", NA, thermo.depth))

# Plot to look at results
ggplot(thermo_depth,mapping=aes(x=datetime,y=-thermo.depth)) +
  geom_line()+
  xlim(as.POSIXct("2020-04-01"), as.POSIXct("2026-04-30"))

## Calculate N2
n2_freq <-rLakeAnalyzer::ts.buoyancy.freq(fcr_temp,seasonal=TRUE)

# Remove NA rows in the key column before joining
n2_freq <- n2_freq |> filter(!is.na(datetime))

#join the N2 and thermocline datasets by datetime
fcr_results_la <- full_join(thermo_depth,n2_freq,by="datetime")

#plot check 
ggplot()+
  geom_line(data=fcr_results_la, aes(x=datetime,y=thermo.depth))

#Calculate schmidt stability
# Create sample bathymetry data
date = seq(as.POSIXct("2020-05-01 00:10:00"), as.POSIXct("2026-04-30 23:59:00"), by = "30 mins")

# Use the bathymetry data on EDI for this. You will also need to cite it. 
fcr_bathy <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/1254/1/f7fa2a06e1229ee75ea39eb586577184?key=yltMpS4UEIk12AvB9L7OL5uRiG0")|>
  filter(Reservoir =="FCR")|>
  dplyr::rename("areas"= SA_m2,
                "depths" = Depth_m)|>
  select(depths, areas)

#calculate schmidt stability using fcr_temp and fcr_bathy
sch_stab <- rLakeAnalyzer::ts.schmidt.stability(fcr_temp,fcr_bathy, na.rm = TRUE)

# Remove NA rows in the key column before joining
sch_stab <- sch_stab |> filter(!is.na(datetime))

# Join data frames together by datetime
fcr_results_la <- full_join(thermo_depth,n2_freq,by="datetime")
fcr_results_la_sch_stab <- full_join(fcr_results_la,sch_stab,by="datetime")

#daily mean
fcr_results_la_sch_stab <- fcr_results_la_sch_stab |>
  mutate(doy = lubridate::yday(datetime),
         Year= lubridate::year(datetime)) |>
  na.omit() |>
  group_by(doy,Year) |>
  summarise(mean_sch = mean(schmidt.stability),
            mean_buo = mean(n2),
            mean_thermodepth = mean(thermo.depth))

#plot thermocline depth
plot <- ggplot()+
  geom_line(data = fcr_results_la_sch_stab,aes(x=doy,y=mean_thermodepth,color=as.factor(Year)),linewidth=1)
plot+scale_y_reverse()

#plot schmidt stability
ggplot()+
  geom_line(data = fcr_results_la_sch_stab,aes(x=doy,y=mean_sch,color=as.factor(Year)),linewidth=1)
############


############
# Let's prepare a dataset containing thermocline depth, buoyancy frequency and schmidt stability using DOY and Year as grouping variables.

la_daily <- fcr_results_la_sch_stab |>
  group_by(Year, doy) |>
  summarise(thermo.depth = mean(mean_thermodepth,na.rm=TRUE),
            n2 = mean(mean_buo,na.rm=TRUE),
            sch_stab = mean(mean_sch,na.rm=TRUE))

# get the date from day of year
la_daily$date <- as.Date(la_daily$doy - 1, origin = paste0(la_daily$Year, "-01-01"))
############



##########
#FINALLY MERGE ALL DATASETS TO CREATE A COMPREHENSIVE DAILY DATASET
#Select the columns of interest to avoid unnecessary columns
#we want dissolved CH4 and CO2 with depths 
ghg_selected <- ghg_fluxes_FCR |>
  select(date,Depth_m, mean_CO2,mean_CH4) |>
  pivot_longer(cols = c(mean_CO2, mean_CH4), names_to = "Category", values_to = "Value") |> 
  pivot_wider(names_from = c(Category,Depth_m), values_from = Value)


#we want single daily value of thermocline depth, buoyancy freq and schmidt stability
la_daily_selected <- la_daily |>
  select(date, thermo.depth, n2, sch_stab)


# ── Rebuild comp_data using met_daily as the complete date spine ───────────────
comp_data <- met_daily |>                          # ← spine: 0 missing days
  left_join(la_daily_selected, by = "date") |>
  left_join(water_temp_9m,     by = "date") |>
  left_join(ec_daily,          by = "date") |>
  left_join(ghg_selected,      by = "date")

comp_data$date <- as.Date(comp_data$date)

# ── Verify ────────────────────────────────────────────────────────────────────
missing_2022_after <- seq(as.Date("2022-01-01"),
                          as.Date("2022-12-31"), by = "day") |>
  (\(x) x[!x %in% as.Date(comp_data$date)])()

cat("Missing 2022 dates after fix:", length(missing_2022_after), "\n")
cat("Total comp_data rows:        ", nrow(comp_data), "\n")

# Check all years
comp_data |>
  mutate(year = year(date)) |>
  group_by(year) |>
  summarise(
    n_days          = n(),
    n_missing_temp5 = sum(is.na(temp_5m)),
    n_missing_ec    = sum(is.na(daily_CH4)),
    n_missing_la    = sum(is.na(DOsat_9m)),
    .groups         = "drop"
  ) |>
  print()

write.csv(comp_data,
          "C:/Users/bibek/Desktop/Bibek_files/5_yrs_FCR_flux_analysis/comp_data.csv",
          row.names = FALSE)

  









