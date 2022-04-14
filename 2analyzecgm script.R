

analysecgm <- function(input,
                         output,
                         filename = "CGM Analysis",
                         aboveexcursionlength = 35,
                         belowexcursionlength = 10,
                         magedef = "1sd",
                         daystart = 6,
                         dayend = 24,
                         format = "rows") {
  
  # Read in tidied data from the input directory. 
  files <- base::list.files(path = input, full.names = TRUE)
  cgmfiles <- 
    base::as.data.frame(base::matrix(nrow = 0,ncol = base::length(files)))
  base::colnames(cgmfiles) <- base::rep("Record",base::length(files))
  # Define the order in which lubridate parses dates.  
  dateparseorder <- c("mdy HM","mdy HMS","mdY HM","mdY HMS","dmy HM","dmy HMS",
                      "dmY HM","dmY HMS","Ymd HM","Ymd HMS","ymd HM","ymd HMS",
                      "Ydm HM","Ydm HMS","ydm HM","ydm HMS")
  allhours <- 0:23
  # read through the input directory, for each file, calculate CGM variables
  
  for (f in 1:base::length(files)) {    
    # Basic variables
    table <- utils::read.csv(files[f],stringsAsFactors = FALSE,na.strings = c("NA",""))
    table$subjectid <- table$subjectid[1]
        table <- unique(table)
    ##obtain studyid, week, and cgmperiod from each file
    cgmfiles["studyid",f] <- table$studyid[1]
    cgmfiles["week",f] <- table$week[1]
    cgmfiles["cgmperiod",f] <- table$cgmperiod[1]
    
    # format table to correct vector (POSIXct, numeric, character et cetera) 
    table$timestamp <- 
      base::as.POSIXct(lubridate::parse_date_time(table$timestamp,
                                                  dateparseorder,tz = "Australia/Sydney"))
    table$sglucose <- base::as.numeric(table$sglucose)
    interval <- as.numeric(900) ## 900 seconds(15 minutes) is the interval period between sensor glucose readings
    cgmfiles["cgm_start_date", f] <- 
      base::as.character(min(table$timestamp,na.rm = T)) ## date of first CGM sensor glucose reading
    
    totaltime <- 
      base::as.numeric(base::difftime(base::max(table$timestamp, na.rm = T),
                                      base::min(table$timestamp,na.rm = T),
                                      units = "secs")) ##difference between first and last sensor glucose reading in seconds
    
    ## calculate percentage of time the CGM is active (there should be a reading every 900 seconds)
    cgmfiles["percent_active",f] <- 
      base::round(((base::length(which(!is.na(table$sglucose)))/(totaltime/interval))*100),2)
    
    table <- table[,-c(1)]
    table <- table[stats::complete.cases(table),] ##only return values which are complete i.e. keeps only the complete rows
   
    ## create a column with glucose readings in mg/dL, to aid some calculations later on
    table$glucosemg <- (table$sglucose) * 18
    
    cgmfiles["num_days_worn",f] <- 
      base::round(base::length(table$sglucose)/(86400/interval), digits = 1)
    cgmfiles["n_total_readings",f] <- 
      base::as.numeric(base::length(base::which(!is.na(table$sglucose))))
    cgmfiles["mean_glucose",f] <- 
      base::round(base::mean(table$sglucose[base::which(!is.na(table$sglucose))],na.rm = T), digits = 2)
    cgmfiles["sd",f] <- 
      base::round(stats::sd(table$sglucose[base::which(!is.na(table$sglucose))]), digits = 2)
    cgmfiles["cv",f] <- 100*
      base::round((stats::sd(table$sglucose[base::which(!is.na(table$sglucose))]))/
      base::mean(table$sglucose[base::which(!is.na(table$sglucose))]), digits = 3)
    cgmfiles["e_a1c",f] <- 
      base::round((2.59 + ((base::mean(table$sglucose[
        base::which(!is.na(table$sglucose))])))) / 1.59,digits = 1)
    cgmfiles["GMI",f] <- base::round(12.71 + (4.70587 * (base::mean(table$sglucose[
      base::which(!is.na(table$sglucose))]))), digits = 2)
    cgmfiles["q1_sensor",f] <- 
      base::as.numeric(base::summary(table$sglucose[
        base::which(!is.na(table$sglucose))])[2])
    cgmfiles["median_sensor",f] <- 
      base::as.numeric(base::summary(table$sglucose[
        base::which(!is.na(table$sglucose))])[3])
    cgmfiles["q3_sensor",f] <- 
      base::as.numeric(base::summary(table$sglucose[
        base::which(!is.na(table$sglucose))])[5])
    cgmfiles["min_glucose",f] <- 
      base::min(table$sglucose[base::which(!is.na(table$sglucose))])
    cgmfiles["max_glucose",f] <- 
      base::max(table$sglucose[base::which(!is.na(table$sglucose))])
    

## Calculate time above range (TAR), time in range (TIR), and time below range (TBR) as per guidelines (2019)
   
    # TBR < 3 mmol/L
    TBR3 <- 
      base::as.numeric(table$sglucose[base::which(!is.na(
        table$sglucose))],length = 1)
    TBR3[TBR3 <= 3] <- 1
    TBR3[TBR3 > 3] <- 0
    BG3.rle <- base::rle(TBR3)
    excursions3 <- 
      base::as.numeric(BG3.rle$lengths[base::which(BG3.rle$values == 1)])
    
    cgmfiles["TBR_3",f] <- 
      base::length(base::which(excursions3 > 
                                 ((belowexcursionlength * 60)/interval)))
    cgmfiles["TBR min 3",f] <- base::sum(TBR3) * (interval/60)
    cgmfiles["TBR percent 3",f] <- 
      ((base::sum(TBR3) * (interval/60))/
         (base::length(table$sglucose) * (interval/60))) * 100
    
    
    # TBR 3.0 - 3.8 mmol/L
    TBR3_38 <- 
      base::as.numeric(table$sglucose[base::which(!is.na(table$sglucose))],
                       length = 1)
    TBR3_38 <- ifelse(TBR3_38 >=3.0 & TBR3_38 <=3.8, 1,0)
    cgmfiles["TBR min 3_3.8",f] <- base::sum(TBR3_38) * (interval/60)
    cgmfiles["TBR percent 3_3.8",f] <- 
      ((base::sum(TBR3_38) * (interval/60))/
         (base::length(table$sglucose) * (interval/60))) * 100   
    
    
    # TIR 3.9 - 10.0 mmol/L
    TIR3_10 <- 
      base::as.numeric(table$sglucose[base::which(!is.na(table$sglucose))],
                       length = 1)
    TIR3_10 <- ifelse(TIR3_10 >=3.9 & TIR3_10 <=10, 1,0)
    cgmfiles["TIR min 3.9_10",f] <- base::sum(TIR3_10) * (interval/60)
    cgmfiles["TIR percent 3.9_10",f] <- 
      ((base::sum(TIR3_10) * (interval/60))/
         (base::length(table$sglucose) * (interval/60))) * 100
    
    
    # TAR 10.1 - 13.9 mmol/L
    TAR10_13 <- 
      base::as.numeric(table$sglucose[base::which(!is.na(table$sglucose))],
                       length = 1)
    TAR10_13 <- ifelse(TAR10_13 >=10.1 & TAR10_13 <=13.9, 1,0)
    cgmfiles["TAR min 10.1_13.9",f] <- base::sum(TAR10_13) * (interval/60)
    cgmfiles["TAR percent 10.1_13.9",f] <- 
      ((base::sum(TAR10_13) * (interval/60))/
         (base::length(table$sglucose) * (interval/60))) * 100
    
     # TAR > 13.9 mmol/L
    TAR13.9 <- 
      base::as.numeric(table$sglucose[base::which(!is.na(
        table$sglucose))],length = 1)
    TAR13.9[TAR13.9 < 13.9] <- 0
    TAR13.9[TAR13.9 >= 13.9] <- 1
    BG13.9.rle <- base::rle(TAR13.9)
    excursions13.9 <- 
      base::as.numeric(BG13.9.rle$lengths[base::which(BG13.9.rle$values == 1)])
    
    cgmfiles["TAR_13.9",f] <- 
      base::length(base::which(excursions13.9 > 
                                 ((aboveexcursionlength * 60)/interval)))
    cgmfiles["TAR min 13.9",f] <- base::sum(TAR13.9) * (interval/60)
    cgmfiles["TAR percent 13.9",f] <- 
      ((base::sum(TAR13.9) * (interval/60))/
         (base::length(table$sglucose) * (interval/60))) * 100
    
    # TAR >= 7.8mmol/L (not in concensus targets, but assoicatied with an increase in macrovascular risk)
    TAR7.8 <- 
      base::as.numeric(table$sglucose[base::which(!is.na(
        table$sglucose))],length = 1)
    TAR7.8[TAR7.8 < 7.8] <- 0
    TAR7.8[TAR7.8 >= 7.8] <- 1
    BG7.8.rle <- base::rle(TAR7.8)
    excursions7.8 <- 
      base::as.numeric(BG7.8.rle$lengths[base::which(BG7.8.rle$values == 1)])
    
    cgmfiles["TAR_7.8",f] <- 
      base::length(base::which(excursions7.8 > 
                                 ((aboveexcursionlength * 60)/interval)))
    cgmfiles["TAR min 7.8",f] <- base::sum(TAR7.8) * (interval/60)
    cgmfiles["TAR percent 7.8",f] <- 
      ((base::sum(TAR7.8) * (interval/60))/
         (base::length(table$sglucose) * (interval/60))) * 100

  ## Find all the "daytime readings"
    daytime_indexes <- 
      base::which(base::as.numeric(base::format(table$timestamp,"%H")) %in% 
                    daystart:dayend)
    daytime_sensor <- table$sglucose[daytime_indexes]
    xaxis <- 
      base::seq(from = 0, length.out = base::length(daytime_sensor),by = 
                  (interval / 60))
  
    
    # Remove NAs if they are present.
    xaxis[base::which(is.na(daytime_sensor))] <- NA
    xaxis <- xaxis[!is.na(xaxis)]
    daytime_sensor <- daytime_sensor[!is.na(daytime_sensor)]
    #
    #calculate AUC using trapezoidal integration for "daytime" readings
    aucs <- pracma::cumtrapz(xaxis,daytime_sensor) ##pracma package
    aucs2 <- MESS::auc(xaxis, daytime_sensor, type = "linear") ## MESS package
    cgmfiles["daytime_auc",f] <- aucs[base::length(daytime_sensor)]
    cgmfiles["daytime_aucMESS",f] <- aucs2
    
    # TIR variables for daytime
    TBR <- ifelse(daytime_sensor < 3, 1,0)
    cgmfiles["day TBR 3 (min)",f] <- base::sum(TBR,na.rm = T) * (interval/60)
    cgmfiles["day TBR 3 (percent)",f] <- 
      (base::sum(TBR,na.rm = T) * (interval/60))/(base::length(daytime_sensor) * (interval/60)) * 100
    
    TBR <- ifelse(daytime_sensor >=3.0 & daytime_sensor <=3.8, 1,0)
    cgmfiles["day TBR 3-3.8 (min)",f] <- base::sum(TBR,na.rm = T) * (interval/60)
    cgmfiles["day TBR 3-3.8 (percent)",f] <- 
      (base::sum(TBR,na.rm = T) * (interval/60))/(base::length(daytime_sensor) * (interval/60)) * 100
    
    
    TIR <- ifelse(daytime_sensor >=3.9 & daytime_sensor <=10, 1,0)
    cgmfiles["day TIR 3.9_10 (min)",f] <- base::sum(TIR,na.rm = T) * (interval/60)
    cgmfiles["day TIR 3.9_10 (percent)",f] <- 
      (base::sum(TIR,na.rm = T) * (interval/60))/(base::length(daytime_sensor) * (interval/60)) * 100
    
  
    TAR <- ifelse(daytime_sensor >=10.1 & daytime_sensor <=13.9, 1,0)
    cgmfiles["day TAR 10.1-13.9 (min)",f] <- base::sum(TAR,na.rm = T) * (interval/60)
    cgmfiles["day TAR 10.1-13.9 (percent)",f] <- 
      (base::sum(TAR,na.rm = T) * (interval/60))/(base::length(daytime_sensor) * (interval/60)) * 100
    
    TAR <- ifelse(daytime_sensor > 13.9, 1,0)
    cgmfiles["day TAR 13.9 (min)",f] <- base::sum(TAR,na.rm = T) * (interval/60)
    cgmfiles["day TAR 13.9 (percnt)",f] <- 
      (base::sum(TAR,na.rm = T) * (interval/60))/(base::length(daytime_sensor) * (interval/60)) * 100
    
    # Other daytime sensor glucose variables.
    cgmfiles["day mean glucose",f] <- 
      base::mean(stats::na.omit(daytime_sensor))
    cgmfiles["day min glucose",f] <- base::min(daytime_sensor)
    cgmfiles["day max glucose",f] <- base::max(daytime_sensor)
    cgmfiles["day sd",f] <- stats::sd(daytime_sensor)
    

    # Nighttime AUC. ? consensus guideline define nocturnal as from midnight to 6am
    nighttime_indexes <- 
      base::which(base::as.numeric(base::format(table$timestamp,"%H")) %in% 
                    allhours[base::which(!(0:23 %in% daystart:dayend))])
    if (length(nighttime_indexes) > 0) {
      nighttime_sensor <- table$sglucose[nighttime_indexes]
      xaxis <- 
        base::seq(from = 0, length.out = base::length(nighttime_indexes),by = 
                    (interval / 60))
      
      # Day/night ratio.
      cgmfiles["day_night_ratio",f] <- 
        base::round(base::length(daytime_sensor)/base::length(nighttime_sensor),1)
      
      # Remove NAs if they are present.
      xaxis[base::which(is.na(nighttime_sensor))] <- NA
      xaxis <- xaxis[!is.na(xaxis)]
      nighttime_sensor <- nighttime_sensor[!is.na(nighttime_sensor)]
      aucs <- pracma::cumtrapz(xaxis,nighttime_sensor)
      cgmfiles["nighttime_auc",f] <- aucs[base::length(nighttime_sensor)]
      
      # TBR, TIR, TAR variables for nighttime
      TBR <- ifelse(nighttime_sensor < 3, 1,0)
      cgmfiles["night TBR 3 (min)",f] <- base::sum(TBR,na.rm = T) * (interval/60)
      cgmfiles["night TBR 3 (percent)",f] <- 
        (base::sum(TBR,na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * (interval/60)) * 100
      
      TBR <- ifelse(nighttime_sensor >=3.0 & nighttime_sensor <=3.8, 1,0)
      cgmfiles["night TBR 3-3.8 (min)",f] <- base::sum(TBR,na.rm = T) * (interval/60)
      cgmfiles["night TBR 3-3.8 (percent)",f] <- 
        (base::sum(TBR,na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * (interval/60)) * 100
      
      
      TIR <- ifelse(nighttime_sensor >=3.9 & nighttime_sensor <=10, 1,0)
      cgmfiles["night TIR 3.9_10 (min)",f] <- base::sum(TIR,na.rm = T) * (interval/60)
      cgmfiles["night TIR 3.9_10 (percent)",f] <- 
        (base::sum(TIR,na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * (interval/60)) * 100
      
      
      TAR <- ifelse(nighttime_sensor >=10.1 & nighttime_sensor <=13.9, 1,0)
      cgmfiles["night TAR 10.1-13.9 (min)",f] <- base::sum(TAR,na.rm = T) * (interval/60)
      cgmfiles["night TAR 10.1-13.9 (percent)",f] <- 
        (base::sum(TAR,na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * (interval/60)) * 100
      
      TAR <- ifelse(nighttime_sensor > 13.9, 1,0)
      cgmfiles["night TAR 13.9 (min)",f] <- base::sum(TAR,na.rm = T) * (interval/60)
      cgmfiles["night TAR 13.9 (percnt)",f] <- 
        (base::sum(TAR,na.rm = T) * (interval/60))/(base::length(nighttime_sensor) * (interval/60)) * 100
      
      # Other nighttime sensor glucose variables.
      cgmfiles["night mean glucose",f] <- 
        base::mean(stats::na.omit(nighttime_sensor))
      cgmfiles["night min glucose",f] <- base::min(nighttime_sensor)
      cgmfiles["night max glucose",f] <- base::max(nighttime_sensor)
      cgmfiles["night sd",f] <- stats::sd(nighttime_sensor) 
      }
    
    
    # Calculate the total AUC (total glucose exposure)
    sensorBG <- base::as.numeric(table$sglucose,length = 1)
    xaxis <- 
      base::seq(from = 0, length.out = base::length(sensorBG),by = 
                  (interval / 60))
    
    # Remove NA in both the xaxis and glucose readings, if there are any
    xaxis[base::which(is.na(sensorBG))] <- NA
    xaxis <- xaxis[!is.na(xaxis)]
    sensorBG <- sensorBG[!is.na(sensorBG)]
    aucs <- pracma::cumtrapz(xaxis,sensorBG)
    cgmfiles["total_auc",f] <- aucs[base::length(sensorBG)]
    
    cgmfiles["average_auc_per_day",f] <-
      base::as.numeric(cgmfiles["total_auc",f]) /
      base::as.numeric(cgmfiles["num_days_worn",f])
    
    # AUC over 7.8.
    sensorover7_8 <- table$sglucose
    sensorover7_8 <- sensorover7_8[sensorover7_8 > 7.8]
    sensorover7_8 <- sensorover7_8[!is.na(sensorover7_8)]
    xaxis <- 
      base::seq(from = 0, length.out = base::length(sensorover7_8),by = 
                  (interval / 60))
    
    # Calculate cumulative AUC, and subtract rectangle where length = 7.8 &
    # width = minutes.
    if (base::length(sensorover7_8) > 1) {
      aucs <- pracma::cumtrapz(xaxis,sensorover7_8)
      aucs <- 
        (aucs[base::length(sensorover7_8)]) - (xaxis[base::length(xaxis)] * 7.8)
      cgmfiles["auc_over_7.8",f] <- aucs
    } else {
      cgmfiles["auc_over_7.8",f] <- 0
    }
    cgmfiles["average_auc_7.8",f] <-
      base::as.numeric(cgmfiles["auc_over_7.8",f]) /
      base::as.numeric(cgmfiles["num_days_worn",f])
    
    
    # AUC over 10.
    sensorover10 <- table$sglucose
    sensorover10 <- sensorover10[sensorover10 >= 10]
    sensorover10 <- sensorover10[!is.na(sensorover10)]
    xaxis <- 
      base::seq(from = 0, length.out = base::length(sensorover10),by = 
                  (interval / 60))
    
    # Calculate cumulative AUC, and subtract recatangle where length = 10 &
    # width = minutes.
    if (base::length(sensorover10) > 1) {
      aucs <- pracma::cumtrapz(xaxis,sensorover10)
      aucs <- 
        (aucs[base::length(sensorover10)]) - (xaxis[base::length(xaxis)] * 10)
      cgmfiles["auc_over_10",f] <- aucs
    } else {
      cgmfiles["auc_over_10",f] <- 0
    }
    cgmfiles["average_auc_10",f] <-
      base::as.numeric(cgmfiles["auc_over_10",f]) /
      base::as.numeric(cgmfiles["num_days_worn",f])
    
    # Calculate MAGE.
    # Smooth data using an exponentially weighted moving average, calculate SD of 
    # unsmoothed data.  
    table$smoothed <- 
      base::as.numeric(zoo::rollapply(zoo::zoo(table$glucosemg), 9, 
                                      function(x) c(1,2,4,8,16,8,4,2,1) %*% 
                                        (x / 46),fill = NA))
    table$smoothed[1:4] <- base::mean(stats::na.omit(table$glucosemg[1:4]))
    table$smoothed[(base::length(table$smoothed)-3):
                     base::length(table$smoothed)] <- 
      base::mean(table$sglucose[(base::length(table$glucosemg)-3):
                                       base::length(table$glucosemg)])
    
    sd <- stats::sd(table$glucosemg)
    # Identify turning points, peaks, and nadirs.
    tpoints <- pastecs::turnpoints(table$smoothed)
    peaks <- base::which(tpoints[["peaks"]] == TRUE)
    pits <- base::which(tpoints[["pits"]] == TRUE)
    # Calculate the difference between each nadir and its following peak. If the     
    # data starts on a peak, remove it. Otherwise remove the final pit to create an 
    # even number of pits and peaks.
    if (tpoints[["firstispeak"]] == TRUE && base::length(peaks) != 
        base::length(pits)) {
      peaks <- peaks[2:base::length(peaks)]
    } else if (tpoints[["firstispeak"]] == FALSE && 
               base::length(peaks) != base::length(pits)) {
      pits <- pits[1:(base::length(pits)-1)]
    }
    differences <- table$glucosemg[peaks] - table$glucosemg[pits]
    
    # Calculate the average of the differences greater than the entire dataset 
    # SD, 2SD, etc.
    if (magedef == "1sd") {
      cgmfiles["r_mage",f] <- 
        (base::mean(stats::na.omit(differences[base::which(differences > sd)]))/18)
    }  ## divide by 18 to convert back to mmol/L
    
    #J-index
    cgmfiles["j_index",f] <- 
      0.324 * (base::mean(table$sglucose, na.rm = T) + 
                 stats::sd(table$sglucose, na.rm = T))^2
    
    # CONGA - continuous overlapping net glycaemic action
    n <- 3600 ## default value for CONGA1 is 60 minutes
    conga.times <- table$timestamp + n
    conga.times <- conga.times[!is.na(conga.times)]
    conga.times <- conga.times[base::order(conga.times)]
    conga.times <- conga.times[base::which(conga.times %in% table$timestamp)]
    begin.times <- conga.times - n
    suppressWarnings(congas <- table$sglucose[base::which(table$timestamp %in% conga.times)] - 
                       table$sglucose[base::which(table$timestamp %in% begin.times)])
    cgmfiles[base::paste0("conga_1"),f] <- stats::sd(congas,na.rm = T)
    
    # MODD.
    table$time <- lubridate::round_date(table$timestamp,"5 minutes")
    table$time <- base::strftime(table$time, format = "%H:%M",tz = "Australia/Sydney")
    moddtable <- 
      base::data.frame(base::matrix(ncol = 2,nrow = 
                                      base::length(unique(table$time))))
    base::colnames(moddtable) <- c("time","mean_differences")
    moddtable$time <- base::unique(table$time)
    
    # For each time, calculate differences (absolute values) and average them.   
    for (r in 1:nrow(moddtable)) {
      moddtable$mean_differences[r] <- 
        base::mean(base::abs(base::diff(table$sglucose[
          base::which(table$time == moddtable$time[r])])))
    }
    # Average the averages.
    cgmfiles["modd",f] <- 
      base::mean(stats::na.omit(moddtable$mean_differences))
    
    # LBGI and HBGI (Kovatchev, Boris P., et al. "Symmetrization of the blood glucose measurement scale and its applications." Diabetes Care 20.11 (1997): 1655-1658)
    a <- 1.026
    b <- 1.861
    y <- 1.794

    table$gluctransform <- y * ((base::log(table$sglucose)^a)-b)
    table$rBG <- 10 * (table$gluctransform^2)
    rl <- table$rBG[base::which(table$gluctransform < 0)]
    rh <- table$rBG[base::which(table$gluctransform > 0)]
    cgmfiles["lbgi",f] <- base::mean(stats::na.omit(rl))
    cgmfiles["hbgi",f] <- base::mean(stats::na.omit(rh))
  }
  
  # Write file.
  cgmfiles <- 
    base::cbind("Variable / Field Name" = rownames(cgmfiles),cgmfiles)
  if (format == "rows") {
    cgmfiles <- base::as.data.frame(base::t(cgmfiles))
    cgmfiles <- cgmfiles[-1,]
  }
  filename <- base::paste(output,"/",filename,".csv",sep = "")
  utils::write.csv(cgmfiles, file = filename,row.names = FALSE)
}
analysecgm("Tidied", "Report")
