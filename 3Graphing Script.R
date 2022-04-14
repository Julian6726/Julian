library(dplyr)
library(tidyverse)

graphcgm <- function(input,
                      output = tempdir(), metadata,
                      tz = "Australia/Sydney",
                      yaxis = c(0,20)){
  
  # Get list of all CGM data files 
  files <- base::list.files(path = input,full.names = TRUE)
  
  # Create a data frame to combine indiviudal CGM values into one dataf rame
  aggregatecgm <- base::data.frame(matrix(ncol = 6,nrow = 0))
  colnames(aggregatecgm) <- c("subjectid","timestamp","sglucose", "cgmperiod", "week", "studyid")
  
  # Run through the input directory, combine all data.   
  for (f in 1:length(files)) {
    cgmdata <- 
      utils::read.csv(files[f],stringsAsFactors = FALSE,header = TRUE,skipNul = TRUE)
    id <- cgmdata$studyid[1]
    cgmdata$subjectid <- id
    aggregatecgm <- rbind(cgmdata,aggregatecgm)
  }
  aggregatecgm$sglucose <- as.numeric(aggregatecgm$sglucose) #convert glucose levels to numeric

  ## read in metadata (study arm and diabetic status for stratification)
  ## here, we read in clinical metadata that is used as a graphing variable. The metadata included columns named study_id, which we use later to join with the aggregatecgm dataframe
  metadata <- read.csv("metadata/metadata.csv", header = TRUE, stringsAsFactors = TRUE)
  metadata[,c(5:6)] <- lapply(metadata[,c(5:6)], factor) #convert columns from integer to factors
  metadata <- metadata[,c(1,5,6)] #keep only the necessary columns
  metadata <- metadata %>% tidyr::drop_na() #remove patients that were not randomised to a study arm (i.e screened but not enroled)
  metadata$study_id <- as.numeric(metadata$study_id) #convert study id to numeric
  metadata$study_id <- sprintf("%02d", as.numeric(metadata$study_id)) # add a leading "0" to study ids
  metadata$study_id <- paste0("DIG", metadata$study_id) # add DIG prefix to study id
  aggregatecgm <- left_join(aggregatecgm, metadata, by = c("studyid" = "study_id"))
  

  
  # Remove rows which contain missing values
  aggregatecgm <- aggregatecgm[stats::complete.cases(aggregatecgm),]
  
  write.csv(aggregatecgm, "aggregatecgm.csv")
  
  # Remove dates, so aggregated CGM data can be ordered by time of day    
  dateparseorder <- c("mdy HM","mdy HMS","mdY HM","mdY HMS","dmy HM","dmy HMS",
                      "dmY HM","dmY HMS","Ymd HM","Ymd HMS","ymd HM","ymd HMS",
                      "Ydm HM","Ydm HMS","ydm HM","ydm HMS")
  aggregatecgm$timestamp <- 
    as.POSIXct(lubridate::parse_date_time(aggregatecgm$timestamp,
                                          dateparseorder,tz = "Australia/Sydney")) ## convert timestamp data to time format
  aggregatecgm$hour <- lubridate::round_date(aggregatecgm$timestamp,"hour") ## round time to the nearest hour
  aggregatecgm$time <- 
    as.POSIXct(strftime(aggregatecgm$timestamp,format = "%H:%M", tz = "Australia/Sydney"),
               format = "%H:%M") #remove date specifier, keep just the time
  aggregatecgm$hourmin <- 
    lubridate::round_date(aggregatecgm$timestamp,"10 minutes") #remove date specifier, round the time to the nearest 10 minute mark
  aggregatecgm$hourmin <- 
    as.POSIXct(strftime(aggregatecgm$hourmin,format = "%H:%M", tz = "Australia/Sydney"),
               format = "%H:%M") # convert to time format
  
  lims <- as.POSIXct(strptime(c("00:00","24:00"),format = "%H:%M")) ## set x-axis limits for graphs
  target <- c(3.9, 10) ##target range for glucose measurements
  range<- as.POSIXct(strptime(c("00:00"),format = "%H:%M")) ## lovation of "target range"


## Graphing using ggplot
  # New facet label names for supp variable
  week.labs <- c("Week 3", "Week 4", "Week 7", "Week 8")
  names(week.labs) <- c("1", "2", "3", "4")
 
    ## Individual average CGM, faceted by week
  for(i in unique(aggregatecgm$studyid)){
    gi<- ggplot2::ggplot(aes(x = time, y = sglucose), data = subset(aggregatecgm, studyid == i))+
      ggplot2::geom_smooth(aes(y = sglucose), method="loess", colour = "blue", se = FALSE)+
      ggplot2::geom_point(ggplot2::aes(y = sglucose),shape =".")+
      ggplot2::geom_ribbon(aes(ymin=3.9, ymax=10, fill="Target Range"), fill="skyblue", alpha=0.5, show.legend = TRUE)+
      ggplot2::labs(title = "Daily Overlay Per Subject (loess Smoothing)", subtitle = (as.character(i)))+
      ggplot2::ylab("Sensor Glucose (mmol/L)")+
      ggplot2::xlab("Time (hour)")+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none")+
      ggplot2::theme(axis.text.x = element_text(face="plain", color="black", size=7, angle=45))+
      ggplot2::theme(axis.text.y = element_text(face="plain", color="black", size=7, angle=0))+
      ggplot2::scale_x_datetime(labels = function(x) format(x, format = "%H:%M"), limits = lims, expand=c(0.07,1)) +
      ggplot2::geom_hline(yintercept = c(3.9,10), linetype="dashed")+ ## horizontal lines for target glucose values >3.9 & < 10
      ggplot2::scale_y_continuous(breaks = sort(c(seq(min(aggregatecgm$sglucose), max(aggregatecgm$sglucose), length.out=5),target)))+
      ggplot2::theme_bw()+
      ggplot2::theme(strip.background=element_rect(fill="white"),
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))+
    ggplot2::facet_grid( ~ week, scales="free", labeller = labeller(week=week.labs))+
      ggplot2::annotate("text", x=range, y=6.95, label = "Target Range", color = "black", size =2, angle = 90, vjust=-0.3, fontface="italic")
     ggplot2::ggsave(filename = sprintf('%sbyweek.png', i), plot = gi,
                        path = output,
                        dpi = 300, width = 17, height = 10, units = "cm")
       }
 
  ## Obtain the "median" (50%) and other percentiles shown as if occurring in a single day, and plot as ambulatory glucose profile (AGP)
  cgmperiod.labs <- c("Week 3-4", "Week 7-8")
  names(cgmperiod.labs) <- c("a", "b")

  for(i in unique(aggregatecgm$studyid)){
    quantiles <- aggregatecgm %>% filter(studyid == i) %>%
      group_by(hourmin, cgmperiod) %>%
      summarise(quantiles = quantile(sglucose, c(0.05, 0.25, 0.5, 0.75, 0.95)), q = c("q0.05", "q0.25", "q0.5", "q0.75", "q0.95"))
    quantiles <- quantiles %>% tidyr::spread(key = q, value = quantiles)
    quantilesa <- quantiles %>% filter(cgmperiod == "a")
    quantilesa$perc5 <-  as.numeric(stats::smooth(quantilesa$q0.05,kind = "3R",twiceit = TRUE))
    quantilesa$perc25 <- as.numeric(stats::smooth(quantilesa$q0.25,kind = "3R",twiceit = TRUE))
    quantilesa$perc50 <- as.numeric(stats::smooth(quantilesa$q0.5,kind = "3R",twiceit = TRUE))
    quantilesa$perc75 <- as.numeric(stats::smooth(quantilesa$q0.75,kind = "3R",twiceit = TRUE))
    quantilesa$perc95 <- as.numeric(stats::smooth(quantilesa$q0.95,kind = "3R",twiceit = TRUE))
    
    quantilesb <- quantiles %>% filter(cgmperiod == "b")
    quantilesb$perc5 <-  as.numeric(stats::smooth(quantilesb$q0.05,kind = "3R",twiceit = TRUE))
    quantilesb$perc25 <- as.numeric(stats::smooth(quantilesb$q0.25,kind = "3R",twiceit = TRUE))
    quantilesb$perc50 <- as.numeric(stats::smooth(quantilesb$q0.5,kind = "3R",twiceit = TRUE))
    quantilesb$perc75 <- as.numeric(stats::smooth(quantilesb$q0.75,kind = "3R",twiceit = TRUE))
    quantilesb$perc95 <- as.numeric(stats::smooth(quantilesb$q0.95,kind = "3R",twiceit = TRUE))
    
    quantiles <- rbind(quantilesa, quantilesb)
    q <- ggplot2::ggplot(data = quantiles, ggplot2::aes(x = hourmin))+ 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc5,ymax = quantiles$perc95, fill = "5th to 95th Percentile"), linetype = 2, alpha = 0.5)+
      ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc25,ymax = quantiles$perc75, fill = "Interquartile Range"), linetype = 3, alpha = 0.5)+
      scale_colour_manual(values=c("dodgerblue","dodgerblue4","black"))+
      geom_line(aes(y=perc50, color = "Median"), colour = "Black", show.legend = TRUE)+ geom_hline(yintercept = c(3.9,10), linetype="solid", colour = "red")+
      ggplot2::scale_x_datetime(labels = function(x) format(x, format = "%H:%M"), date_breaks = "3 hours", limits = lims, expand=c(0,0))+
      theme_bw()+
      ggplot2::theme(legend.title = element_blank(), strip.background=element_rect(fill="white"))+
      ggplot2::ylab("Sensor Glucose (mmol/L)")+
      ggplot2::xlab("Time (hour)")+
      ggplot2::labs(title = "Ambulatory Glucose Profile", subtitle = as.character(i))+
      ggplot2::facet_grid(rows = vars(cgmperiod), labeller = labeller(cgmperiod=cgmperiod.labs))
    ggplot2::ggsave(filename = sprintf('%sAGP.png', i), plot = q,
                    path = output,
                    dpi = 300, width = 17, height = 10, units = "cm")
   }

  
## AGP for all patients, stratified by CGM period and DM status
   aggregatecgm0 <- aggregatecgm %>% group_by(hourmin, cgmperiod, post_tx_dm) %>% 
    summarise(quantiles = quantile(sglucose, c(0.05, 0.25, 0.5, 0.75, 0.95)), q = c("q0.05", "q0.25", "q0.5", "q0.75", "q0.95"))
  quantiles <- aggregatecgm0 %>% tidyr::spread(key = q, value = quantiles)
  quantiles0a <- quantiles %>% filter(cgmperiod == "a" & post_tx_dm == "0")
  quantiles0b <- quantiles %>% filter(cgmperiod == "b" & post_tx_dm == "0")
  quantiles1a <- quantiles %>% filter(cgmperiod == "a" & post_tx_dm == "1")
  quantiles1b <- quantiles %>% filter(cgmperiod == "b" & post_tx_dm == "1")
  
  quantiles0a$perc5 <-  as.numeric(stats::smooth(quantiles0a$q0.05,kind = "3R",twiceit = TRUE))
  quantiles0a$perc25 <- as.numeric(stats::smooth(quantiles0a$q0.25,kind = "3R",twiceit = TRUE))
  quantiles0a$perc50 <- as.numeric(stats::smooth(quantiles0a$q0.5,kind = "3R",twiceit = TRUE))
  quantiles0a$perc75 <- as.numeric(stats::smooth(quantiles0a$q0.75,kind = "3R",twiceit = TRUE))
  quantiles0a$perc95 <- as.numeric(stats::smooth(quantiles0a$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles0b$perc5 <-  as.numeric(stats::smooth(quantiles0b$q0.05,kind = "3R",twiceit = TRUE))
  quantiles0b$perc25 <- as.numeric(stats::smooth(quantiles0b$q0.25,kind = "3R",twiceit = TRUE))
  quantiles0b$perc50 <- as.numeric(stats::smooth(quantiles0b$q0.5,kind = "3R",twiceit = TRUE))
  quantiles0b$perc75 <- as.numeric(stats::smooth(quantiles0b$q0.75,kind = "3R",twiceit = TRUE))
  quantiles0b$perc95 <- as.numeric(stats::smooth(quantiles0b$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles1a$perc5 <-  as.numeric(stats::smooth(quantiles1a$q0.05,kind = "3R",twiceit = TRUE))
  quantiles1a$perc25 <- as.numeric(stats::smooth(quantiles1a$q0.25,kind = "3R",twiceit = TRUE))
  quantiles1a$perc50 <- as.numeric(stats::smooth(quantiles1a$q0.5,kind = "3R",twiceit = TRUE))
  quantiles1a$perc75 <- as.numeric(stats::smooth(quantiles1a$q0.75,kind = "3R",twiceit = TRUE))
  quantiles1a$perc95 <- as.numeric(stats::smooth(quantiles1a$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles1b$perc5 <-  as.numeric(stats::smooth(quantiles1b$q0.05,kind = "3R",twiceit = TRUE))
  quantiles1b$perc25 <- as.numeric(stats::smooth(quantiles1b$q0.25,kind = "3R",twiceit = TRUE))
  quantiles1b$perc50 <- as.numeric(stats::smooth(quantiles1b$q0.5,kind = "3R",twiceit = TRUE))
  quantiles1b$perc75 <- as.numeric(stats::smooth(quantiles1b$q0.75,kind = "3R",twiceit = TRUE))
  quantiles1b$perc95 <- as.numeric(stats::smooth(quantiles1b$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles <- rbind(quantiles0a, quantiles0b, quantiles1a, quantiles1b)
  
  ############## 
  
  ##plot and facet by study_arm and cgmperiod
  cgmperiod.labs <- c("Week 3-4", "Week 7-8")
  names(cgmperiod.labs) <- c("a", "b")
  post_tx_dm.labs <- c("Non Diabetic", "Diabetic")
  names(post_tx_dm.labs) <- c("0", "1")
  
  k <- ggplot2::ggplot(data = quantiles, ggplot2::aes(x = hourmin, y=perc50))+ 
    ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc5,ymax = quantiles$perc95, fill = "5th to 95th Percentile"),linetype = 2,alpha = 0.5)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc25,ymax = quantiles$perc75, fill = "Interquartile Range"), linetype = 3, alpha = 0.5)+
    ggplot2::scale_colour_manual(values=c("dodgerblue","dodgerblue4","black"))+
    ggplot2::geom_line(data=quantiles,ggplot2::aes(y=perc50, color = "Median"), colour = "Blue", show.legend = TRUE)+ 
    ggplot2::geom_hline(yintercept = c(3.9,10), linetype="solid", colour = "red")+
    ggplot2::scale_x_datetime(labels = function(x) format(x, format = "%H:%M"), date_breaks = "3 hours", limits = lims, expand=c(0,0))+
    ggplot2::theme_bw()+
    ggplot2::scale_y_continuous(breaks=c(0,3.9,5,7.5,10,12.5,15))+
    ggplot2::theme(legend.title = ggplot2::element_blank(), strip.background=ggplot2::element_rect(fill="white"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size =8),
                   axis.text.y = ggplot2::element_text(size = 8))+
    ggplot2::ylab("Sensor Glucose (mmol/L)")+
    ggplot2::xlab("Time (hour)")+
    ggplot2::labs(title = "Ambulatory Glucose Profile", subtitle = "All Patients")+
    theme(panel.spacing.x = unit(-0.015, "lines"))+
    ggplot2::facet_grid(rows = vars(post_tx_dm), cols=vars(cgmperiod), labeller = labeller(cgmperiod=cgmperiod.labs, post_tx_dm=post_tx_dm.labs))
  ggplot2::ggsave(filename = sprintf('%sAGP.png', "all patients"), plot = k,
                  path = output,
                  dpi = 300, width = 17, height = 10, units = "cm")
  
  
  

  ##plot all aggregate results on the one graph, for non diabetic patients, stratified by study arm and cgm period
  aggregatecgm0 <- aggregatecgm %>% dplyr::filter(post_tx_dm == "0")
  aggregatecgm0 <- aggregatecgm0 %>% group_by(hourmin, study_arm, cgmperiod) %>% 
    summarise(quantiles = quantile(sglucose, c(0.05, 0.25, 0.5, 0.75, 0.95)), q = c("q0.05", "q0.25", "q0.5", "q0.75", "q0.95"))
  quantiles <- aggregatecgm0 %>% tidyr::spread(key = q, value = quantiles)
  quantiles0a <- quantiles %>% filter(cgmperiod == "a" & study_arm == "0")
  quantiles0b <- quantiles %>% filter(cgmperiod == "b" & study_arm == "0")
  quantiles1a <- quantiles %>% filter(cgmperiod == "a" & study_arm == "1")
  quantiles1b <- quantiles %>% filter(cgmperiod == "b" & study_arm == "1")
  
  quantiles0a$perc5 <-  as.numeric(stats::smooth(quantiles0a$q0.05,kind = "3R",twiceit = TRUE))
  quantiles0a$perc25 <- as.numeric(stats::smooth(quantiles0a$q0.25,kind = "3R",twiceit = TRUE))
  quantiles0a$perc50 <- as.numeric(stats::smooth(quantiles0a$q0.5,kind = "3R",twiceit = TRUE))
  quantiles0a$perc75 <- as.numeric(stats::smooth(quantiles0a$q0.75,kind = "3R",twiceit = TRUE))
  quantiles0a$perc95 <- as.numeric(stats::smooth(quantiles0a$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles0b$perc5 <-  as.numeric(stats::smooth(quantiles0b$q0.05,kind = "3R",twiceit = TRUE))
  quantiles0b$perc25 <- as.numeric(stats::smooth(quantiles0b$q0.25,kind = "3R",twiceit = TRUE))
  quantiles0b$perc50 <- as.numeric(stats::smooth(quantiles0b$q0.5,kind = "3R",twiceit = TRUE))
  quantiles0b$perc75 <- as.numeric(stats::smooth(quantiles0b$q0.75,kind = "3R",twiceit = TRUE))
  quantiles0b$perc95 <- as.numeric(stats::smooth(quantiles0b$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles1a$perc5 <-  as.numeric(stats::smooth(quantiles1a$q0.05,kind = "3R",twiceit = TRUE))
  quantiles1a$perc25 <- as.numeric(stats::smooth(quantiles1a$q0.25,kind = "3R",twiceit = TRUE))
  quantiles1a$perc50 <- as.numeric(stats::smooth(quantiles1a$q0.5,kind = "3R",twiceit = TRUE))
  quantiles1a$perc75 <- as.numeric(stats::smooth(quantiles1a$q0.75,kind = "3R",twiceit = TRUE))
  quantiles1a$perc95 <- as.numeric(stats::smooth(quantiles1a$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles1b$perc5 <-  as.numeric(stats::smooth(quantiles1b$q0.05,kind = "3R",twiceit = TRUE))
  quantiles1b$perc25 <- as.numeric(stats::smooth(quantiles1b$q0.25,kind = "3R",twiceit = TRUE))
  quantiles1b$perc50 <- as.numeric(stats::smooth(quantiles1b$q0.5,kind = "3R",twiceit = TRUE))
  quantiles1b$perc75 <- as.numeric(stats::smooth(quantiles1b$q0.75,kind = "3R",twiceit = TRUE))
  quantiles1b$perc95 <- as.numeric(stats::smooth(quantiles1b$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles <- rbind(quantiles0a, quantiles0b, quantiles1a, quantiles1b)
  
  ############## 
  
  ##plot and facet by study_arm and cgmperiod
  cgmperiod.labs <- c("Week 3-4", "Week 7-8")
  names(cgmperiod.labs) <- c("a", "b")
  study_arm.labs <- c("Standard Care", "Inulin")
  names(study_arm.labs) <- c("0", "1")

  h <- ggplot2::ggplot(data = quantiles, ggplot2::aes(x = hourmin, y=perc50))+ 
    ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc5,ymax = quantiles$perc95, fill = "5th to 95th Percentile"),linetype = 2,alpha = 0.5)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc25,ymax = quantiles$perc75, fill = "Interquartile Range"),linetype = 3,  alpha = 0.5)+
    ggplot2::scale_colour_manual(values=c("dodgerblue","dodgerblue4","black"))+
    ggplot2::geom_line(data=quantiles,ggplot2::aes(y=perc50, color = "Median"), colour = "Blue", show.legend = TRUE)+ 
    ggplot2::geom_hline(yintercept = c(3.9,10), linetype="solid", colour = "red")+
    ggplot2::scale_x_datetime(labels = function(x) format(x, format = "%H:%M"), date_breaks = "3 hours", limits = lims, expand=c(0,0))+
    ggplot2::theme_bw()+
    ggplot2::scale_y_continuous(breaks=c(0,3.9,5,7.5,10,12.5,15))+
    ggplot2::theme(legend.title = ggplot2::element_blank(), strip.background=ggplot2::element_rect(fill="white"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size =8),
                   axis.text.y = ggplot2::element_text(size = 8))+
    ggplot2::ylab("Sensor Glucose (mmol/L)")+
    ggplot2::xlab("Time (hour)")+
    ggplot2::labs(title = "Ambulatory Glucose Profile", subtitle = "Non-diabetics")+
    theme(panel.spacing.x = unit(-0.015, "lines"))+
    ggplot2::facet_grid(rows = vars(study_arm), cols=vars(cgmperiod), labeller = labeller(cgmperiod=cgmperiod.labs, study_arm=study_arm.labs))
  ggplot2::ggsave(filename = sprintf('%sAGP.png', "non-diabetic by study arm"), plot = h,
                  path = output,
                  dpi = 300, width = 17, height = 10, units = "cm")
  
  
  ##plot all aggregate results on the one graph, for non diabetic patients, and cgm period only
  aggregatecgm0 <- aggregatecgm %>% dplyr::filter(post_tx_dm == "0")
  aggregatecgm0 <- aggregatecgm0 %>% group_by(hourmin, cgmperiod) %>% 
    summarise(quantiles = quantile(sglucose, c(0.05, 0.25, 0.5, 0.75, 0.95)), q = c("q0.05", "q0.25", "q0.5", "q0.75", "q0.95"))
  quantiles <- aggregatecgm0 %>% tidyr::spread(key = q, value = quantiles)
  quantiles0a <- quantiles %>% filter(cgmperiod == "a")
  quantiles0b <- quantiles %>% filter(cgmperiod == "b")
  quantiles1a <- quantiles %>% filter(cgmperiod == "a")
  quantiles1b <- quantiles %>% filter(cgmperiod == "b")
  
  quantiles0a$perc5 <-  as.numeric(stats::smooth(quantiles0a$q0.05,kind = "3R",twiceit = TRUE))
  quantiles0a$perc25 <- as.numeric(stats::smooth(quantiles0a$q0.25,kind = "3R",twiceit = TRUE))
  quantiles0a$perc50 <- as.numeric(stats::smooth(quantiles0a$q0.5,kind = "3R",twiceit = TRUE))
  quantiles0a$perc75 <- as.numeric(stats::smooth(quantiles0a$q0.75,kind = "3R",twiceit = TRUE))
  quantiles0a$perc95 <- as.numeric(stats::smooth(quantiles0a$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles0b$perc5 <-  as.numeric(stats::smooth(quantiles0b$q0.05,kind = "3R",twiceit = TRUE))
  quantiles0b$perc25 <- as.numeric(stats::smooth(quantiles0b$q0.25,kind = "3R",twiceit = TRUE))
  quantiles0b$perc50 <- as.numeric(stats::smooth(quantiles0b$q0.5,kind = "3R",twiceit = TRUE))
  quantiles0b$perc75 <- as.numeric(stats::smooth(quantiles0b$q0.75,kind = "3R",twiceit = TRUE))
  quantiles0b$perc95 <- as.numeric(stats::smooth(quantiles0b$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles1a$perc5 <-  as.numeric(stats::smooth(quantiles1a$q0.05,kind = "3R",twiceit = TRUE))
  quantiles1a$perc25 <- as.numeric(stats::smooth(quantiles1a$q0.25,kind = "3R",twiceit = TRUE))
  quantiles1a$perc50 <- as.numeric(stats::smooth(quantiles1a$q0.5,kind = "3R",twiceit = TRUE))
  quantiles1a$perc75 <- as.numeric(stats::smooth(quantiles1a$q0.75,kind = "3R",twiceit = TRUE))
  quantiles1a$perc95 <- as.numeric(stats::smooth(quantiles1a$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles1b$perc5 <-  as.numeric(stats::smooth(quantiles1b$q0.05,kind = "3R",twiceit = TRUE))
  quantiles1b$perc25 <- as.numeric(stats::smooth(quantiles1b$q0.25,kind = "3R",twiceit = TRUE))
  quantiles1b$perc50 <- as.numeric(stats::smooth(quantiles1b$q0.5,kind = "3R",twiceit = TRUE))
  quantiles1b$perc75 <- as.numeric(stats::smooth(quantiles1b$q0.75,kind = "3R",twiceit = TRUE))
  quantiles1b$perc95 <- as.numeric(stats::smooth(quantiles1b$q0.95,kind = "3R",twiceit = TRUE))
  
  quantiles <- rbind(quantiles0a, quantiles0b, quantiles1a, quantiles1b)
  
  cgmperiod.labs <- c("Week 3-4", "Week 7-8")
  names(cgmperiod.labs) <- c("a", "b")
  study_arm.labs <- c("Standard Care", "Inulin")
  names(study_arm.labs) <- c("0", "1")
  
  hh <- ggplot2::ggplot(data = quantiles, ggplot2::aes(x = hourmin, y=perc50))+ 
    ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc5,ymax = quantiles$perc95, fill = "5th to 95th Percentile"),linetype = 2,alpha = 0.5)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc25,ymax = quantiles$perc75, fill = "Interquartile Range"),linetype = 3,  alpha = 0.5)+
    ggplot2::scale_colour_manual(values=c("dodgerblue","dodgerblue4","black"))+
    ggplot2::geom_line(data=quantiles,ggplot2::aes(y=perc50, color = "Median"), colour = "Blue", show.legend = TRUE)+ 
    ggplot2::geom_hline(yintercept = c(3.9,10), linetype="solid", colour = "red")+
    ggplot2::scale_x_datetime(labels = function(x) format(x, format = "%H:%M"), date_breaks = "3 hours", limits = lims, expand=c(0,0))+
    ggplot2::theme_bw()+
    ggplot2::scale_y_continuous(breaks=c(0,3.9,5,7.5,10,12.5,15))+
    ggplot2::theme(legend.title = ggplot2::element_blank(), strip.background=ggplot2::element_rect(fill="white"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size =8),
                   axis.text.y = ggplot2::element_text(size = 8))+
    ggplot2::ylab("Sensor Glucose (mmol/L)")+
    ggplot2::xlab("Time (hour)")+
    ggplot2::labs(title = "Ambulatory Glucose Profile", subtitle = "Non-diabetic")+
    theme(panel.spacing.x = unit(-0.015, "lines"))+
    ggplot2::facet_grid(rows=vars(cgmperiod), labeller = labeller(cgmperiod=cgmperiod.labs, study_arm=study_arm.labs))
  ggplot2::ggsave(filename = sprintf('%sAGP.png', "non-diabetic"), plot = hh,
                  path = output,
                  dpi = 300, width = 17, height = 10, units = "cm")
  

    ##plot all aggregate results on the one graph, for DIABETIC patients, stratified by study arm and cgm period
    aggregatecgm1 <- aggregatecgm %>% dplyr::filter(post_tx_dm == "1")
    aggregatecgm1 <- aggregatecgm1 %>% group_by(hourmin, study_arm, cgmperiod) %>% 
      summarise(quantiles = quantile(sglucose, c(0.05, 0.25, 0.5, 0.75, 0.95)), q = c("q0.05", "q0.25", "q0.5", "q0.75", "q0.95"))
    quantiles <- aggregatecgm1 %>% tidyr::spread(key = q, value = quantiles)
    quantiles0a <- quantiles %>% filter(cgmperiod == "a" & study_arm == "0")
    quantiles0b <- quantiles %>% filter(cgmperiod == "b" & study_arm == "0")
    quantiles1a <- quantiles %>% filter(cgmperiod == "a" & study_arm == "1")
    quantiles1b <- quantiles %>% filter(cgmperiod == "b" & study_arm == "1")
    
    quantiles0a$perc5 <-  as.numeric(stats::smooth(quantiles0a$q0.05,kind = "3R",twiceit = TRUE))
    quantiles0a$perc25 <- as.numeric(stats::smooth(quantiles0a$q0.25,kind = "3R",twiceit = TRUE))
    quantiles0a$perc50 <- as.numeric(stats::smooth(quantiles0a$q0.5,kind = "3R",twiceit = TRUE))
    quantiles0a$perc75 <- as.numeric(stats::smooth(quantiles0a$q0.75,kind = "3R",twiceit = TRUE))
    quantiles0a$perc95 <- as.numeric(stats::smooth(quantiles0a$q0.95,kind = "3R",twiceit = TRUE))
    
    quantiles0b$perc5 <-  as.numeric(stats::smooth(quantiles0b$q0.05,kind = "3R",twiceit = TRUE))
    quantiles0b$perc25 <- as.numeric(stats::smooth(quantiles0b$q0.25,kind = "3R",twiceit = TRUE))
    quantiles0b$perc50 <- as.numeric(stats::smooth(quantiles0b$q0.5,kind = "3R",twiceit = TRUE))
    quantiles0b$perc75 <- as.numeric(stats::smooth(quantiles0b$q0.75,kind = "3R",twiceit = TRUE))
    quantiles0b$perc95 <- as.numeric(stats::smooth(quantiles0b$q0.95,kind = "3R",twiceit = TRUE))
    
    quantiles1a$perc5 <-  as.numeric(stats::smooth(quantiles1a$q0.05,kind = "3R",twiceit = TRUE))
    quantiles1a$perc25 <- as.numeric(stats::smooth(quantiles1a$q0.25,kind = "3R",twiceit = TRUE))
    quantiles1a$perc50 <- as.numeric(stats::smooth(quantiles1a$q0.5,kind = "3R",twiceit = TRUE))
    quantiles1a$perc75 <- as.numeric(stats::smooth(quantiles1a$q0.75,kind = "3R",twiceit = TRUE))
    quantiles1a$perc95 <- as.numeric(stats::smooth(quantiles1a$q0.95,kind = "3R",twiceit = TRUE))
    
    quantiles1b$perc5 <-  as.numeric(stats::smooth(quantiles1b$q0.05,kind = "3R",twiceit = TRUE))
    quantiles1b$perc25 <- as.numeric(stats::smooth(quantiles1b$q0.25,kind = "3R",twiceit = TRUE))
    quantiles1b$perc50 <- as.numeric(stats::smooth(quantiles1b$q0.5,kind = "3R",twiceit = TRUE))
    quantiles1b$perc75 <- as.numeric(stats::smooth(quantiles1b$q0.75,kind = "3R",twiceit = TRUE))
    quantiles1b$perc95 <- as.numeric(stats::smooth(quantiles1b$q0.95,kind = "3R",twiceit = TRUE))
    
    quantiles <- rbind(quantiles0a, quantiles0b, quantiles1a, quantiles1b)
    
    ##plot and facet by study_arm and cgmperiod
    cgmperiod.labs <- c("Week 3-4", "Week 7-8")
    names(cgmperiod.labs) <- c("a", "b")
    study_arm.labs <- c("Standard Care", "Inulin")
    names(study_arm.labs) <- c("0", "1")
    j <- ggplot2::ggplot(data = quantiles, ggplot2::aes(x = hourmin, y=perc50))+ 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc5,ymax = quantiles$perc95, fill = "5th to 95th Percentile"),linetype = 2,alpha = 0.5)+
      ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc25,ymax = quantiles$perc75, fill = "Interquartile Range"), linetype = 3, alpha = 0.5)+
      ggplot2::geom_line(ggplot2::aes(y=perc50, color = perc50), colour = "blue", show.legend = TRUE)+ 
      ggplot2::geom_hline(yintercept = c(3.9,10), linetype="solid", colour = "red")+
      ggplot2::scale_y_continuous(breaks=c(0,3.9,7.5,10,15,20))+
      ggplot2::scale_x_datetime(labels = function(x) format(x, format = "%H:%M"), date_breaks = "3 hours", limits = lims, expand=c(0,0))+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.title = ggplot2::element_blank(), strip.background=ggplot2::element_rect(fill="white"),
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8), 
                     axis.text.y = ggplot2::element_text(size = 8),legend.position = "right")+
      ggplot2::ylab("Sensor Glucose (mmol/L)")+
      ggplot2::xlab("Time (hour)")+
      ggplot2::labs(title = "Ambulatory Glucose Profile", subtitle = "Diabetics")+
      ggplot2::facet_grid(rows = vars(study_arm), cols=vars(cgmperiod), 
                          labeller = labeller(cgmperiod=cgmperiod.labs, study_arm=study_arm.labs))+
      theme(panel.spacing.x = unit(-0.015, "lines"))
      ggplot2::ggsave(filename = sprintf('%sAGP.png', "diabetic by study arm"), plot = j,
                    path = output,
                    dpi = 300, width = 17, height = 10, units = "cm")  
      
      ##plot all aggregate results on the one graph, for DIABETIC patients, stratify by cgm period 
      aggregatecgm0 <- aggregatecgm %>% dplyr::filter(post_tx_dm == "1")
      aggregatecgm0 <- aggregatecgm0 %>% group_by(hourmin, cgmperiod) %>% 
        summarise(quantiles = quantile(sglucose, c(0.05, 0.25, 0.5, 0.75, 0.95)), q = c("q0.05", "q0.25", "q0.5", "q0.75", "q0.95"))
      quantiles <- aggregatecgm0 %>% tidyr::spread(key = q, value = quantiles)
      quantiles0a <- quantiles %>% filter(cgmperiod == "a")
      quantiles0b <- quantiles %>% filter(cgmperiod == "b")
      quantiles1a <- quantiles %>% filter(cgmperiod == "a")
      quantiles1b <- quantiles %>% filter(cgmperiod == "b")
      
      quantiles0a$perc5 <-  as.numeric(stats::smooth(quantiles0a$q0.05,kind = "3R",twiceit = TRUE))
      quantiles0a$perc25 <- as.numeric(stats::smooth(quantiles0a$q0.25,kind = "3R",twiceit = TRUE))
      quantiles0a$perc50 <- as.numeric(stats::smooth(quantiles0a$q0.5,kind = "3R",twiceit = TRUE))
      quantiles0a$perc75 <- as.numeric(stats::smooth(quantiles0a$q0.75,kind = "3R",twiceit = TRUE))
      quantiles0a$perc95 <- as.numeric(stats::smooth(quantiles0a$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles0b$perc5 <-  as.numeric(stats::smooth(quantiles0b$q0.05,kind = "3R",twiceit = TRUE))
      quantiles0b$perc25 <- as.numeric(stats::smooth(quantiles0b$q0.25,kind = "3R",twiceit = TRUE))
      quantiles0b$perc50 <- as.numeric(stats::smooth(quantiles0b$q0.5,kind = "3R",twiceit = TRUE))
      quantiles0b$perc75 <- as.numeric(stats::smooth(quantiles0b$q0.75,kind = "3R",twiceit = TRUE))
      quantiles0b$perc95 <- as.numeric(stats::smooth(quantiles0b$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles1a$perc5 <-  as.numeric(stats::smooth(quantiles1a$q0.05,kind = "3R",twiceit = TRUE))
      quantiles1a$perc25 <- as.numeric(stats::smooth(quantiles1a$q0.25,kind = "3R",twiceit = TRUE))
      quantiles1a$perc50 <- as.numeric(stats::smooth(quantiles1a$q0.5,kind = "3R",twiceit = TRUE))
      quantiles1a$perc75 <- as.numeric(stats::smooth(quantiles1a$q0.75,kind = "3R",twiceit = TRUE))
      quantiles1a$perc95 <- as.numeric(stats::smooth(quantiles1a$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles1b$perc5 <-  as.numeric(stats::smooth(quantiles1b$q0.05,kind = "3R",twiceit = TRUE))
      quantiles1b$perc25 <- as.numeric(stats::smooth(quantiles1b$q0.25,kind = "3R",twiceit = TRUE))
      quantiles1b$perc50 <- as.numeric(stats::smooth(quantiles1b$q0.5,kind = "3R",twiceit = TRUE))
      quantiles1b$perc75 <- as.numeric(stats::smooth(quantiles1b$q0.75,kind = "3R",twiceit = TRUE))
      quantiles1b$perc95 <- as.numeric(stats::smooth(quantiles1b$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles <- rbind(quantiles0a, quantiles0b, quantiles1a, quantiles1b)
      
      cgmperiod.labs <- c("Week 3-4", "Week 7-8")
      names(cgmperiod.labs) <- c("a", "b")
      study_arm.labs <- c("Standard Care", "Inulin")
      names(study_arm.labs) <- c("0", "1")
      
      jj <- ggplot2::ggplot(data = quantiles, ggplot2::aes(x = hourmin, y=perc50))+ 
        ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc5,ymax = quantiles$perc95, fill = "5th to 95th Percentile"),linetype = 2,alpha = 0.5)+
        ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc25,ymax = quantiles$perc75, fill = "Interquartile Range"),linetype = 3,  alpha = 0.5)+
        ggplot2::scale_colour_manual(values=c("dodgerblue","dodgerblue4","black"))+
        ggplot2::geom_line(data=quantiles,ggplot2::aes(y=perc50, color = "Median"), colour = "Blue", show.legend = TRUE)+ 
        ggplot2::geom_hline(yintercept = c(3.9,10), linetype="solid", colour = "red")+
        ggplot2::scale_x_datetime(labels = function(x) format(x, format = "%H:%M"), date_breaks = "3 hours", limits = lims, expand=c(0,0))+
        ggplot2::theme_bw()+
        ggplot2::scale_y_continuous(breaks=c(0,3.9,5,7.5,10,12.5,15))+
        ggplot2::theme(legend.title = ggplot2::element_blank(), strip.background=ggplot2::element_rect(fill="white"),
                       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size =8),
                       axis.text.y = ggplot2::element_text(size = 8))+
        ggplot2::ylab("Sensor Glucose (mmol/L)")+
        ggplot2::xlab("Time (hour)")+
        ggplot2::labs(title = "Ambulatory Glucose Profile", subtitle = "Diabetic")+
        theme(panel.spacing.x = unit(-0.015, "lines"))+
        ggplot2::facet_grid(rows=vars(cgmperiod), labeller = labeller(cgmperiod=cgmperiod.labs, study_arm=study_arm.labs))
      ggplot2::ggsave(filename = sprintf('%sAGP.png', "diabetic"), plot = jj,
                      path = output,
                      dpi = 300, width = 17, height = 10, units = "cm")
      
      ##plot all aggregate results on the one graph, for non diabetic patients, and cgm period only
      aggregatecgm0 <- aggregatecgm %>% dplyr::filter(post_tx_dm == "0")
      aggregatecgm0 <- aggregatecgm0 %>% group_by(hourmin, cgmperiod) %>% 
        summarise(quantiles = quantile(sglucose, c(0.05, 0.25, 0.5, 0.75, 0.95)), q = c("q0.05", "q0.25", "q0.5", "q0.75", "q0.95"))
      quantiles <- aggregatecgm0 %>% tidyr::spread(key = q, value = quantiles)
      quantiles0a <- quantiles %>% filter(cgmperiod == "a")
      quantiles0b <- quantiles %>% filter(cgmperiod == "b")
      quantiles1a <- quantiles %>% filter(cgmperiod == "a")
      quantiles1b <- quantiles %>% filter(cgmperiod == "b")
      
      quantiles0a$perc5 <-  as.numeric(stats::smooth(quantiles0a$q0.05,kind = "3R",twiceit = TRUE))
      quantiles0a$perc25 <- as.numeric(stats::smooth(quantiles0a$q0.25,kind = "3R",twiceit = TRUE))
      quantiles0a$perc50 <- as.numeric(stats::smooth(quantiles0a$q0.5,kind = "3R",twiceit = TRUE))
      quantiles0a$perc75 <- as.numeric(stats::smooth(quantiles0a$q0.75,kind = "3R",twiceit = TRUE))
      quantiles0a$perc95 <- as.numeric(stats::smooth(quantiles0a$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles0b$perc5 <-  as.numeric(stats::smooth(quantiles0b$q0.05,kind = "3R",twiceit = TRUE))
      quantiles0b$perc25 <- as.numeric(stats::smooth(quantiles0b$q0.25,kind = "3R",twiceit = TRUE))
      quantiles0b$perc50 <- as.numeric(stats::smooth(quantiles0b$q0.5,kind = "3R",twiceit = TRUE))
      quantiles0b$perc75 <- as.numeric(stats::smooth(quantiles0b$q0.75,kind = "3R",twiceit = TRUE))
      quantiles0b$perc95 <- as.numeric(stats::smooth(quantiles0b$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles1a$perc5 <-  as.numeric(stats::smooth(quantiles1a$q0.05,kind = "3R",twiceit = TRUE))
      quantiles1a$perc25 <- as.numeric(stats::smooth(quantiles1a$q0.25,kind = "3R",twiceit = TRUE))
      quantiles1a$perc50 <- as.numeric(stats::smooth(quantiles1a$q0.5,kind = "3R",twiceit = TRUE))
      quantiles1a$perc75 <- as.numeric(stats::smooth(quantiles1a$q0.75,kind = "3R",twiceit = TRUE))
      quantiles1a$perc95 <- as.numeric(stats::smooth(quantiles1a$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles1b$perc5 <-  as.numeric(stats::smooth(quantiles1b$q0.05,kind = "3R",twiceit = TRUE))
      quantiles1b$perc25 <- as.numeric(stats::smooth(quantiles1b$q0.25,kind = "3R",twiceit = TRUE))
      quantiles1b$perc50 <- as.numeric(stats::smooth(quantiles1b$q0.5,kind = "3R",twiceit = TRUE))
      quantiles1b$perc75 <- as.numeric(stats::smooth(quantiles1b$q0.75,kind = "3R",twiceit = TRUE))
      quantiles1b$perc95 <- as.numeric(stats::smooth(quantiles1b$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles <- rbind(quantiles0a, quantiles0b, quantiles1a, quantiles1b)
      
      hh <- ggplot2::ggplot(data = quantiles, ggplot2::aes(x = hourmin, y=perc50, color = cgmperiod)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc5,ymax = quantiles$perc95, fill = "5th to 95th Percentile"), linetype = 2, alpha = 0.2) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc25,ymax = quantiles$perc75, fill = "Interquartile Range"), linetype = 3, alpha = 0.4) +
        ggplot2::geom_line(data=quantiles,ggplot2::aes(y=perc50, color = cgmperiod), show.legend = TRUE)+ 
        ggplot2::scale_color_brewer(palette = "Set1", labels = c("Weeks 3-4", "Weeks 7-8"))+
        ggplot2::geom_hline(yintercept = c(3.9,10), linetype="solid", colour = "red")+
        ggplot2::scale_x_datetime(labels = function(x) format(x, format = "%H:%M"), date_breaks = "3 hours", limits = lims, expand=c(0,0))+
        ggplot2::theme_bw()+
        ggplot2::scale_y_continuous(breaks=c(0,3.9,5,7.5,10,12.5,15))+
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size =8),
                       axis.text.y = ggplot2::element_text(size = 8), legend.position = "top")+
        ggplot2::ylab("Sensor Glucose (mmol/L)")+
        ggplot2::xlab("Time (hour)")+
        ggplot2::labs(title = "Ambulatory Glucose Profile", subtitle = "Non-diabetic Patients")+
        theme(panel.spacing.x = unit(-0.015, "lines"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))
      ggplot2::ggsave(filename = sprintf('%sAGP.png', "non-diabetic overlay"), plot = hh,
                      path = output,
                      dpi = 300, width = 17, height = 10, units = "cm")

      ##plot all aggregate results on the one graph, for diabetic patients, and cgm period only
      aggregatecgm0 <- aggregatecgm %>% dplyr::filter(post_tx_dm == "1")
      aggregatecgm0 <- aggregatecgm0 %>% group_by(hourmin, cgmperiod) %>% 
        summarise(quantiles = quantile(sglucose, c(0.05, 0.25, 0.5, 0.75, 0.95)), q = c("q0.05", "q0.25", "q0.5", "q0.75", "q0.95"))
      quantiles <- aggregatecgm0 %>% tidyr::spread(key = q, value = quantiles)
      quantiles0a <- quantiles %>% filter(cgmperiod == "a")
      quantiles0b <- quantiles %>% filter(cgmperiod == "b")
      quantiles1a <- quantiles %>% filter(cgmperiod == "a")
      quantiles1b <- quantiles %>% filter(cgmperiod == "b")
      
      quantiles0a$perc5 <-  as.numeric(stats::smooth(quantiles0a$q0.05,kind = "3R",twiceit = TRUE))
      quantiles0a$perc25 <- as.numeric(stats::smooth(quantiles0a$q0.25,kind = "3R",twiceit = TRUE))
      quantiles0a$perc50 <- as.numeric(stats::smooth(quantiles0a$q0.5,kind = "3R",twiceit = TRUE))
      quantiles0a$perc75 <- as.numeric(stats::smooth(quantiles0a$q0.75,kind = "3R",twiceit = TRUE))
      quantiles0a$perc95 <- as.numeric(stats::smooth(quantiles0a$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles0b$perc5 <-  as.numeric(stats::smooth(quantiles0b$q0.05,kind = "3R",twiceit = TRUE))
      quantiles0b$perc25 <- as.numeric(stats::smooth(quantiles0b$q0.25,kind = "3R",twiceit = TRUE))
      quantiles0b$perc50 <- as.numeric(stats::smooth(quantiles0b$q0.5,kind = "3R",twiceit = TRUE))
      quantiles0b$perc75 <- as.numeric(stats::smooth(quantiles0b$q0.75,kind = "3R",twiceit = TRUE))
      quantiles0b$perc95 <- as.numeric(stats::smooth(quantiles0b$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles1a$perc5 <-  as.numeric(stats::smooth(quantiles1a$q0.05,kind = "3R",twiceit = TRUE))
      quantiles1a$perc25 <- as.numeric(stats::smooth(quantiles1a$q0.25,kind = "3R",twiceit = TRUE))
      quantiles1a$perc50 <- as.numeric(stats::smooth(quantiles1a$q0.5,kind = "3R",twiceit = TRUE))
      quantiles1a$perc75 <- as.numeric(stats::smooth(quantiles1a$q0.75,kind = "3R",twiceit = TRUE))
      quantiles1a$perc95 <- as.numeric(stats::smooth(quantiles1a$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles1b$perc5 <-  as.numeric(stats::smooth(quantiles1b$q0.05,kind = "3R",twiceit = TRUE))
      quantiles1b$perc25 <- as.numeric(stats::smooth(quantiles1b$q0.25,kind = "3R",twiceit = TRUE))
      quantiles1b$perc50 <- as.numeric(stats::smooth(quantiles1b$q0.5,kind = "3R",twiceit = TRUE))
      quantiles1b$perc75 <- as.numeric(stats::smooth(quantiles1b$q0.75,kind = "3R",twiceit = TRUE))
      quantiles1b$perc95 <- as.numeric(stats::smooth(quantiles1b$q0.95,kind = "3R",twiceit = TRUE))
      
      quantiles <- rbind(quantiles0a, quantiles0b, quantiles1a, quantiles1b)
      
      kk <- ggplot2::ggplot(data = quantiles, ggplot2::aes(x = hourmin, y=perc50, color = cgmperiod)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc5,ymax = quantiles$perc95, fill = "5th to 95th Percentile"), linetype = 2, alpha = 0.2) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = quantiles$perc25,ymax = quantiles$perc75, fill = "Interquartile Range"), linetype = 3, alpha = 0.4) +
        ggplot2::geom_line(data=quantiles,ggplot2::aes(y=perc50, color = cgmperiod), show.legend = TRUE)+ 
        ggplot2::scale_color_brewer(palette = "Set1", labels = c("Weeks 3-4", "Weeks 7-8"))+
        ggplot2::geom_hline(yintercept = c(3.9,10), linetype="solid", colour = "red")+
        ggplot2::scale_x_datetime(labels = function(x) format(x, format = "%H:%M"), date_breaks = "3 hours", limits = lims, expand=c(0,0))+
        ggplot2::theme_bw()+
        ggplot2::scale_y_continuous(breaks=c(0,3.9,5,7.5,10,12.5,15,20,25))+
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size =8),
                       axis.text.y = ggplot2::element_text(size = 8), legend.position = "top")+
        ggplot2::ylab("Sensor Glucose (mmol/L)")+
        ggplot2::xlab("Time (hour)")+
        ggplot2::labs(title = "Ambulatory Glucose Profile", subtitle = "Diabetic Patients")+
        theme(panel.spacing.x = unit(-0.015, "lines"))+
        guides(color=guide_legend(override.aes=list(fill=NA)))
      ggplot2::ggsave(filename = sprintf('%sAGP.png', "diabetic overlay"), plot = kk,
                      path = output,
                      dpi = 300, width = 17, height = 10, units = "cm")
 }
 
graphcgm("Practice","Graphs","metadata")






