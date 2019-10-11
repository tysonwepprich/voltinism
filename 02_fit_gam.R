
source('01_data_prep.R')


# # for parallel model fitting
# library(foreach) # for parallelized loops
# library(doParallel)
# 
# # run a species or subset for testing
# # allspecies <- allspecies[c(81:84), ]
# 
# ncores <- 8
# if(ncores > (parallel::detectCores() / 2)){
#   ncores <- parallel::detectCores() / 2
# }
# cl <- makePSOCKcluster(ncores)
# registerDoParallel(cl)
# 
# mcoptions <- list(preschedule = FALSE)

# # foreach loop
# outfiles <- foreach(sp = 1:nrow(allspecies),
#                     .combine='c',
#                     .packages= c("mgcv", "dplyr", "tidyr", "purrr",
#                                  "lubridate"),
#                     .export = c("data", "surveys", "covdata", "allspecies", "gdd"),
#                     .inorder = FALSE,
#                     .options.multicore = mcoptions) %dopar% {

# 42 attempts at mixture models, not all will work
mvspecies <- read.csv("data/speciestraits.csv", header = TRUE) %>% 
  filter(mv_analysis == 1) %>% 
  droplevels.data.frame()

sp <- 32

species <- mvspecies$CommonName[sp]
maxvolt <- mvspecies$Voltinism[sp]


# TODO: bootstrap from the beginning? Could do it by unique survey or SiteYear


counts <- data %>% 
  filter(CommonName == species) %>% 
  mutate(DOY = yday(SiteDate),
         Year = year(SiteDate))

#get unique surveys, including those where species not counted (but were at one time)
survs <- surveys %>% 
  # filter(year(SiteDate) %in% unique(counts$Year)) %>% 
  filter(SiteID %in% unique(counts$SiteID)) %>% 
  mutate(Year = year(SiteDate))

#Add zeros to surveys when species not counted during a survey
test <- left_join(survs, counts, by = c("SiteID", "SiteDate", "Week", "SeqID", "Year"))
counts <- test[, c("SiteID", "SiteDate", "Week", "SeqID", "Total", "Year")]
counts$Total <- plyr::mapvalues(counts$Total, from = NA, to = 0)
counts <- left_join(counts, covdata)
counts$temperature[which(counts$temperature < 50)] <- NA
counts$duration[which(counts$duration == 0)] <- NA

# trying to add GDD instead of ordinal date
counts <- counts %>% 
  left_join(gdd) %>% 
  group_by(SiteID, Year) %>% 
  mutate(SurvPerYear = length(unique(SeqID)),
         YearTotal = sum(Total))

# what if cutoff for inclusion is really open?
# could filter out sites with lower effort later
# fits with UKBMS approach to use all data for GAM for imputation, 
# then filter sites with too many missing surveys in 2nd step

dat <- counts %>% filter(YearTotal >= 1, SurvPerYear >= 10)

# mod <- list()
  
  dat$Year <- as.factor(as.character(dat$Year))
  dat$region <- as.factor(as.character(dat$region))
  dat$SiteID <- as.factor(as.character(dat$SiteID))
  dat$SiteYear <- as.factor(paste(dat$SiteID, dat$Year, sep = "_"))
  dat$RegYear <- as.factor(paste(dat$region, dat$Year, sep = "_"))
  dat <- as.data.frame(dat)
  
  dat <- dat[which(!is.na(dat$AccumDD)), ]
  
  
  
  dd_dist <- rep(dat$AccumDD, dat$Total)
  mm <- Mclust(dd_dist, G = 1:maxvolt, modelNames = "E")
  summary(mm, parameters = TRUE)
  
  # dd_dist <- rep(dat$DOY, dat$Total)
  # mm <- Mclust(dd_dist, G = 1:maxvolt, modelNames = "E")
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
    geom_point(alpha = .2)
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  outlist <- list()
  for (ngen in 1:length(mm$parameters$pro)){
  
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
           SiteYear = as.factor(SiteYear),
           RegYear = as.factor(RegYear))
  
  mod <- bam(GenTotal ~ 
               # s(listlength) + #using listlength depresses model predictions in general
               s(AccumDD, bs = "cr", k = 10, m = 2) +
               s(AccumDD, region, k = 5, bs = "fs", m = 1) +
               s(AccumDD, RegYear, k = 5, bs = "fs", m = 1) +
               # te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(3, 10), d = c(2, 1), by = Year) +
               s(SiteID, bs = "re") +
               s(Year, bs = "re") +
               s(SiteYear, bs = "re"),
             # family = nb(theta = NULL, link = "log"),
             family = poisson(link = "log"),
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
    filter(GenTotal > 0) %>% 
    dplyr::select(AccumDD, DOY, region, GenTotal) %>% 
    filter(complete.cases(.))
  
  outs$region <- factor(outs$region, levels = c("NW", "NE", "SW", "CN"))
  
  
  gamplt <- ggplot(preds, aes(x = AccumDD, y = adjY, group = SiteYear, color = SiteYearGDD)) +
    geom_path(aes(alpha = log(SiteYearTotal))) + 
    geom_point(data = outs, aes(x = AccumDD, y = GenTotal), alpha = .2, inherit.aes = FALSE) +
    scale_color_viridis() + 
    facet_wrap(~region, ncol = 2) +
    geom_rug(data = outs, aes(x = AccumDD), sides="b", inherit.aes = FALSE,  alpha = 0.3) +
    ggtitle("seasonal phenology modeled on degree-day scale") +
    labs(color = "Total degree-days\n for site and year") +
    labs(x = "Degree-days accumulated (5/30C thresholds)") +
    labs(y = "Scaled phenology (model predictions)")
  gamplt
  