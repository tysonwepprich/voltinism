# Preparing data, source this first
# Ohio Lepidopterist volunteers have contributed these observations
# They ask that the data be used for research and that the group be acknowledged

# packages to load 
library(mclust)
library(lubridate)
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(viridis)

theme_set(theme_bw(base_size = 14)) 

data <- readr::read_csv("data/data.trim.csv") %>% 
  mutate(SiteID = formatC(SiteID.x, width = 3, format = "d", flag = "0"),
         SiteDate = lubridate::ymd(SiteDate))

# Two cryptic species are not distinguished in the monitoring
data$CommonName[which(data$CommonName == "Spring/Summer Azure")] <- "Azures"

# 1995 was a pilot year
data <- data %>% 
  filter(year(SiteDate) >= 1996, year(SiteDate) <= 2016)

# Filter out unidentified species
allspecies <- data %>% 
  filter(CommonName %in% unique(CommonName)[1:122]) %>% 
  group_by(CommonName, CombinedLatin) %>% 
  summarise(n = sum(Total)) %>% 
  arrange(n)


surveys <- distinct(data[, c("SeqID", "SiteID", "SiteDate", "Week")])

# Covariates for surveys
# Listlength is # of species observed, often used as catch-all covariate for effort/weather/season/etc.
covdata <- data %>%
  group_by(SeqID) %>%
  summarise(listlength = length(which(unique(CommonName) %in% allspecies$CommonName)),
            temperature = mean(c(StartTemp, EndTemp), na.rm = TRUE),
            duration = duration[1]) %>%
  distinct() %>% 
  left_join(surveys)

# impute missing duration values (due to misaligned time entries, etc.)
# only one still is NA because never reported end time for surveys.
covdata <- covdata %>% 
  mutate(Year = year(SiteDate)) %>% 
  group_by(SiteID, Year) %>% 
  mutate(duration = ifelse(is.na(duration) == TRUE | duration == 0,
                           median(duration, na.rm = TRUE),
                           duration))

# Geographic coordinates for sites, address/transect start (no info in data about route)
sites <- read.csv("data/OHsites2018update.txt") %>% 
  mutate(SiteID = formatC(Name, width = 3, format = "d", flag = "0"))

# Growing degree-day data for each site using Daymet interpolations of daily tmax/tmin
# https://daymet.ornl.gov/getdata
# 
gdd <- readRDS("data/dailyDD.rds")

# 5C/30C thresholds used (instead of more typical 10C for a generic lower threshold)
# if 10C used, some sites/years would have zero accumulation before monitoring started April 1
gdd <- left_join(gdd, sites) %>% 
  dplyr::select(SiteID, SiteDate, degday530, chill0, lat, lon, maxT, minT) %>% 
  mutate(Year = year(SiteDate),
         DOY = yday(SiteDate)) %>% 
  group_by(SiteID, Year) %>% 
  arrange(DOY) %>% 
  mutate(AccumDD = cumsum(degday530),
         AccumChill = cumsum(chill0))

siteGDD <- gdd %>%
  group_by(SiteID, lat, lon) %>% 
  # mutate(season = case_when(month(SiteDate) %in% c(12, 1, 2) ~ "winter",
  #                           month(SiteDate) %in% c(3, 4, 5) ~ "spring",
  #                           month(SiteDate) %in% c(6, 7, 8) ~ "summer",
  #                           month(SiteDate) %in% c(9, 10, 11) ~ "fall"),
  #        )
  filter(DOY == 365) %>%
  summarise(meanGDD = mean(AccumDD),
            meanChill = mean(AccumChill))

# many ways to cluster sites, but using lat/lon is simplest 
# wanted 4 regions for plotting simplicity
sitemod <- densityMclust(siteGDD[,c(2:3)], G = 1:4, modelNames = "EVV")
siteGDD$region <- as.character(sitemod$classification)
# visualize regional clusters
a <- ggplot(data = siteGDD, aes(x = lon, y = lat, group = region, color = region)) + geom_point()
a

siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                  to = c("NE", "NW", "CN", "SW"))
gdd <- gdd %>% 
  left_join(siteGDD[, c("SiteID", "region")])

# graph of temperature trends over time
# summary of gdd/mean temp for manuscript
# Daymet starts at 1980

temperature <- gdd %>%
  filter(Year >= 1980) %>% # or 1980
  group_by(SiteID, Year) %>%
  summarise(meantemp = mean((maxT + minT / 2), na.rm = TRUE),
            meangdd5 = max(AccumDD)) %>%
  ungroup() %>%
  mutate(SiteID = as.factor(as.character(SiteID))) %>%
  group_by(SiteID) %>%
  mutate(meansitetemp = mean(meantemp),
         meansitegdd5 = mean(meangdd5)) %>%
  group_by(Year) %>%
  mutate(meanyeartemp = mean(meantemp),
         meanyeargdd5 = mean(meangdd5))

tm <- lme4::lmer(meantemp ~ Year + (1 | SiteID), data = temperature)
tg <- gam(meantemp ~ s(Year, k = 6) + s(SiteID, bs = "re"), data = temperature)
plot(tg)

summary(tm)

temp <- gdd %>% 
  # filter(Year >= 1996) %>%
  group_by(Year) %>% 
  summarise(meantemp = mean((maxT + minT / 2), na.rm = TRUE),
            meangdd5 = max(AccumDD))

tm <- lm(meangdd5 ~ Year, data = temp)
summary(tm)
