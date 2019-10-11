devtools::install_github("gavinsimpson/gratia")
library(gratia)

dat <- dat %>% filter(year(SiteDate) == 2011) %>% droplevels.data.frame() %>% 
  group_by(Week) %>% 
  mutate(zll = listlength - mean(listlength))
  

mod <- gam(Total ~ s(listlength) +
                     te(lat, lon, DOY, bs = c("tp", "cr"), k = c(8, 40), d = c(2, 1)) +
                     s(SiteID, bs = "re"),
               family = nb(theta = NULL, link = "log"),
               data = dat,
                   method = "REML",
                   control = list(maxit = 500))

mod_nb <- gam(Total ~ s(listlength) + s(DOY, k = 5) +
                te(AccumDD, SiteID, bs = c("cr", "re")),
              family = nb(theta = NULL, link = "log"),
              data = dat,
              method = "REML",
              control = list(maxit = 500))

mod_site <- bam(Total ~s(listlength) +
                  te(AccumDD, bs = "cr", k = 30, by = region) +
                  s(Year, bs = "re") +
                  s(SiteID, bs = "re") +
                  s(SiteYear, DOY, bs = "fs", k = 5, m = 1),
                family = nb(theta = NULL, link = "log"),
                data = dat,
                method = "fREML",
                control = list(maxit = 500))


mod_year <- bam(Total ~s(listlength) +
                  te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1), by = Year) +
                  s(SiteYear, bs = "re") +
                  s(Year, bs = "re") +
                  s(SiteID, bs = "re"),
                family = nb(theta = NULL, link = "log"),
                data = dat,
                method = "fREML",
                control = list(maxit = 500))


datplot <- dat %>% filter(SurvPerYear > 15, YearTotal > 10) %>% droplevels.data.frame()
p <- ggplot(datplot, aes(x = AccumDD, y = Total, group = SiteID)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 20), method.args = list(family = nb(theta = NULL, link = "log"))) +
  facet_wrap(~SiteID, ncol = 5, scales = "free")
p



counts <- dat %>% filter(YearTotal >= 5)

preds <- gdd %>% 
  mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
  filter(SiteYear %in% unique(counts$SiteYear)) %>%
  group_by(SiteID, Year) %>% 
  mutate(SiteYearGDD = max(AccumDD)) %>% 
  filter(DOY %in% seq(90, 305, 1)) %>% 
  ungroup() %>% 
  mutate(RegYear = as.factor(paste(region, Year, sep = "_" )),
         Year = as.factor(as.character(Year)),
         SiteID = as.factor(SiteID),
         SiteYear = as.factor(SiteYear),
         listlength = mean(dat$listlength))

preds$adjY <- predict.gam(object = mod, newdata = preds, type="response")

preds <- preds %>% 
  group_by(SiteYear) %>%
  mutate(Gamma = adjY / as.vector(sum(adjY)),
         SiteYearTotal = sum(adjY)) %>%
  ungroup() %>% 
  filter(adjY > 0) #%>% 
  # mutate(ztotal = ((.7 - .1) * (SiteYearTotal - min(SiteYearTotal))
  #                  /(max(SiteYearTotal) - min(SiteYearTotal))) + .1) %>% 
  # group_by(region) %>% 
  # filter(SiteYear %in% sample(unique(SiteYear), 10, replace = TRUE))

preds$region <- factor(preds$region, levels = c("NW", "NE", "SW", "CN"))

# outliers
outs <- counts %>% 
  ungroup() %>% 
  filter(Total > 0) %>% 
  dplyr::select(AccumDD, DOY, region) %>% 
  filter(complete.cases(.))

outs$region <- factor(outs$region, levels = c("NW", "NE", "SW", "CN"))


gamplt <- ggplot(preds, aes(x = DOY, y = Gamma, group = SiteYear, color = SiteYearGDD)) +
  geom_path(aes(alpha = log(SiteYearTotal))) + 
  scale_color_viridis() + 
  facet_wrap(~region, ncol = 2) +
  geom_rug(data = outs, aes(x = DOY), sides="b", inherit.aes = FALSE,  alpha = 0.3) +
  ggtitle("seasonal phenology modeled on degree-day scale") +
  labs(color = "Total degree-days\n for site and year") +
  labs(x = "Degree-days accumulated (5/30C thresholds)") +
  labs(y = "Scaled phenology (model predictions)")
gamplt

dotplt <- ggplot(counts, aes(x = AccumDD, y = Total, color = SiteID)) +
  geom_point() +
  facet_wrap(~region, ncol = 2)
dotplt



dd_dist <- rep(dat$AccumDD, dat$Total)
mm <- Mclust(dd_dist, G = 1:4, modelNames = "E")
summary(mm, parameters = TRUE)

# dd_dist <- rep(dat$DOY, dat$Total)
# mm <- Mclust(dd_dist, G = 1:5)
# summary(mm, parameters = TRUE)

# assign generations
dat$row <- 1:nrow(dat)
dat2 <- dat[rep(dat$row, dat$Total), ]
dat2 <- rbind(dat2, dat[which(dat$Total == 0), ])

quickclass <- function(x){rmultinom(n = 1, size = 1, prob = x)}

# add DOY to compare mixture models
mmpred <- predict(mm, newdata = dat2$AccumDD)
res <- as.data.frame(mmpred$z)
resclass <- apply(X = res, MARGIN = 1, FUN  = quickclass)
resclass <- data.frame(t(resclass))
dat2$gen <- apply(X = resclass, MARGIN = 1, FUN = function(x) which(x == 1))

dat3 <- dat2 %>% 
  group_by(SiteID, SiteDate, gen) %>% 
  summarise(GenTotal = length(Total)) %>% 
  left_join(dat2) %>% 
  distinct()

dat3$GenTotal[which(dat3$Total == 0)] <- 0

ggplot(dat3, aes(x = AccumDD, y = GenTotal, color = as.factor(gen))) +
  geom_point(alpha = .1)

'%!in%' <- function(x,y)!('%in%'(x,y))
ngen <- 3

dat4 <- dat3 %>% filter(gen == ngen) %>% droplevels.data.frame()
# assign unselected gen counts/zeros as all zero counts for gen of interest?
dat4a <- dat3 %>% 
  filter(gen != ngen) %>% 
  filter(SeqID %!in% dat4$SeqID) %>% 
  mutate(GenTotal = 0) %>% 
  filter(AccumDD < max(dat4$AccumDD) + 200,
         AccumDD > min(dat4$AccumDD) - 200)

dat4 <- bind_rows(dat4,dat4a) %>% 
  ungroup() %>% 
  mutate(SiteID = as.factor(SiteID),
         Year = as.factor(Year),
         SiteYear = as.factor(SiteYear))

mod <- bam(GenTotal ~ 
              s(listlength) +
              # ti(DOY) +
              # ti(DOY, listlength) +
              s(AccumDD, bs = "cr", k = 10, by = Year) +
             # te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1), by = Year) +
              s(SiteID, bs = "re") +
              s(Year, bs = "re") +
             s(SiteYear, bs = "re"),
              # s(AccumDD, SiteYear, bs = "fs"),
            family = nb(theta = NULL, link = "log"),
           # family = poisson(link = "log"),
           data = dat4,
           discrete = TRUE, nthreads = 2,
           method = "fREML",
           control = list(maxit = 500))



counts <- dat4 #%>% filter(YearTotal >= 5)

preds <- gdd %>% 
  mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
  filter(SiteYear %in% unique(counts$SiteYear)) %>%
  group_by(SiteID, Year) %>% 
  mutate(SiteYearGDD = max(AccumDD)) %>% 
  filter(AccumDD >= min(counts$AccumDD[counts$GenTotal > 0])) %>% 
  filter(AccumDD <= max(counts$AccumDD[counts$GenTotal > 0])) %>% 
  # filter(DOY %in% seq(min(counts$DOY), max(counts$DOY), 1)) %>% 
  ungroup() %>% 
  mutate(RegYear = as.factor(paste(region, Year, sep = "_" )),
         Year = as.factor(as.character(Year)),
         SiteID = as.factor(SiteID),
         SiteYear = as.factor(SiteYear),
         listlength = mean(dat4$listlength))

preds$adjY <- predict.gam(object = mod, newdata = preds, type="response")

preds <- preds %>% 
  group_by(SiteYear) %>%
  mutate(Gamma = adjY / as.vector(sum(adjY)),
         SiteYearTotal = sum(adjY)) %>%
  ungroup() %>% 
  filter(adjY > 0) #%>% 
# mutate(ztotal = ((.7 - .1) * (SiteYearTotal - min(SiteYearTotal))
#                  /(max(SiteYearTotal) - min(SiteYearTotal))) + .1) %>% 
# group_by(region) %>% 
# filter(SiteYear %in% sample(unique(SiteYear), 10, replace = TRUE))

preds$region <- factor(preds$region, levels = c("NW", "NE", "SW", "CN"))

# outliers
outs <- counts %>% 
  ungroup() %>% 
  filter(Total > 0) %>% 
  dplyr::select(AccumDD, DOY, region) %>% 
  filter(complete.cases(.))

outs$region <- factor(outs$region, levels = c("NW", "NE", "SW", "CN"))


gamplt <- ggplot(preds, aes(x = AccumDD, y = adjY, group = SiteYear, color = SiteYearGDD)) +
  geom_path(aes(alpha = log(SiteYearTotal))) + 
  scale_color_viridis() + 
  facet_wrap(~region, ncol = 2) +
  geom_rug(data = outs, aes(x = AccumDD), sides="b", inherit.aes = FALSE,  alpha = 0.3) +
  ggtitle("seasonal phenology modeled on degree-day scale") +
  labs(color = "Total degree-days\n for site and year") +
  labs(x = "Degree-days accumulated (5/30C thresholds)") +
  labs(y = "Scaled phenology (model predictions)")
gamplt




# fireflies
ff <- readxl::read_xlsx("../../../Downloads/FFW_Data_2016.xlsx") %>% 
  filter(Longitude < -50, Latitude > 20, Year < 2017) %>% 
  mutate(DOY = yday(`Observation Date`)) %>% 
  mutate(LatGroup = cut(Latitude, breaks = 4))

ggplot(ff, aes(x = Longitude, y = Latitude, color = DOY)) +
  geom_point() +
  scale_color_viridis() +
  facet_wrap(~Year)

ggplot(ff, aes(x = DOY, group = Year, color = Year)) +
  geom_density() +
  facet_wrap(~LatGroup, ncol = 2, scales = "free")

dat <- ff %>% filter(LatGroup == unique(LatGroup)[3])
dd_dist <- rep(dat$DOY, dat$`Number Seen in 10 s`)
mm <- Mclust(dd_dist, G = 1:4, modelNames = "E")
summary(mm, parameters = TRUE)

