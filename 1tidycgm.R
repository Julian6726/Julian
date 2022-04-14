

library(dplyr)
library(tidyverse)
library(lubridate)


## create a function to tidy the raw CGM data (in CSV file), downloaded from LibreView website
tidycgm <- function(input,
                      output,
                      removegaps = TRUE,
                      gapfill = TRUE,
                      maximumgap = 20,
                    week = TRUE) {

  ## The system is directed to the folder where the unformatted .csv files are stored, and to where the new files will be deposited once formatting is complete.
  # directory.
  files <- base::list.files(path = input,full.names = TRUE)
  base::dir.create(output,showWarnings = FALSE)
  dateparseorder <- c("mdy HM","mdy HMS","mdY HM","mdY HMS","dmy HM","dmy HMS",
                      "dmY HM","dmY HMS","Ymd HM","Ymd HMS","ymd HM","ymd HMS",
                      "Ydm HM","Ydm HMS","ydm HM","ydm HMS")
  
  # Read in data from .csv files.  
  for (f in 1:base::length(files)) {
      table <- utils::read.csv(files[f],
                               stringsAsFactors = FALSE,
                               header = TRUE,
                               na.strings = "")

    
          # Format columns to remove extraneous data and blank fields
      id <- table[1,1] ##locate the subjects "id" from the unformatted file
      base::colnames(table) <- table[2,] ## rename columns based on row 2, rather than the header
      table <- table[-c(1:2),] ## removes the first 2 rows of data
      table <- table[,c("Device Timestamp","Historic Glucose mmol/L")] ## keep only the columns named "Device Timestamp","Historic Glucose mmol/L"
      base::colnames(table) <- c('timestamp','sglucose') ## rename the columns to coding friendly variables
  
    
    # If necessary, remove rows with no data.
    if (NA %in% table$timestamp) {
      table <- table[-c(base::which(is.na(table$timestamp))),]
    }
    table <- na.omit(table) ##Remove sensor readings that record NA values
    
    ##format the "timestamp" data into POSIXct format, representing calender dates and times
    table$timestamp <- base::as.POSIXct(lubridate::parse_date_time(table$timestamp,dateparseorder, tz="Australia/Sydney"))
    
    ##format the "sglucose" column to numeric format
    table$sglucose <- base::suppressWarnings(base::as.numeric(table$sglucose))
    ##order the readings by "timestamp" in chronological order
    table <- table[base::order(table$timestamp),]
    
    ## extract the date and time the first glucose measurement is recorded
    recordstart <- 
      base::strftime(table$timestamp[min(which(!is.na(table$sglucose)))],
                     format = "%m/%d/%Y %T")
    ## extract the date and time the last glucose measurement is recorded
    recordstop <- 
      base::strftime(table$timestamp[length(table$timestamp)],
                     format = "%m/%d/%Y %T")
    
    # Set interval based on mode of differenc between timestamp (the most frequent interval)
    interval <- pracma::Mode(base::diff(base::as.numeric(table$timestamp)))
    
    # Clean data (optional).
    if (removegaps == TRUE) {
      # Remove first rows without sensor glucose data.      
      if (is.na(table$sglucose[1])) {
        table <- 
          table[-c(1:base::min(base::which(!is.na(table$sglucose))) - 1),]
      }
      # Remove first 4 hours of data based on timestamp. Add 14,400 seconds (4 hours) 
      # to first timestamp. 
      hour4 <- base::as.numeric(table$timestamp[1]) + 14400
      # Determine which row contains the timestamp closest to hour4, remove all rows 
      # up to and including that row.    
      table <- 
        table[-c(1:(which(abs(as.numeric(table$timestamp) - hour4) == 
                            min(abs(as.numeric(table$timestamp) - hour4)))[1])),]
     
      # Fill in small sensor glucose data gaps using interpolated values
      if (gapfill == TRUE) {
        table$sglucose <- zoo::na.approx(table$sglucose,na.rm = FALSE,
                                              maxgap = (maximumgap*60)/interval)
      }
      # If remaining gaps are larger than the maximum, remove the 24 chunk containing 
      # the gap.  
      repeat(
        if (NA %in% table$sglucose) {
          # Determine the start time for the sensor data gap.          
          startNA <- 
            base::as.numeric(table$timestamp[base::min(base::which(is.na(
              table$sglucose)))])
          # Add 24 hours minus one recording interval.          
          hour24 <- startNA + (86400 - interval)
          table <- table[-c(base::suppressWarnings(base::which(base::abs(
            base::as.numeric(table$timestamp) - startNA) == base::min(base::abs(
              base::as.numeric(table$timestamp) - startNA))):(base::which(
                base::abs(base::as.numeric(
                  table$timestamp) - hour24) == base::min(base::abs(
                    base::as.numeric(table$timestamp) - hour24)))))),]
        } else if (!(NA %in% table$sglucose)) {
          break()
        }
      )
      if (base::length(table$timestamp) == 0) {
        stop(base::paste("File '",files[f],"' does not have enough data and 
                         cannot be processed with the current settings.",
                         sep = ""))
      }
      # Trim end of data so it is in 24 hour chunks.
      seconds <- 
        ((base::as.numeric(base::floor(table$timestamp[base::length(
          table$timestamp)] - table$timestamp[1]))) * 86400) - interval
      table <- 
        table[-c(base::which(table$timestamp > 
                               (table$timestamp[1] + seconds))),]
      if ((1 - base::as.numeric(table$timestamp[base::length(
        table$timestamp)] - table$timestamp[1])%%1) > 0.1) {
        seconds <- ((base::as.numeric(base::floor(table$timestamp[base::length(
          table$timestamp)] - table$timestamp[1]))) * 86400) - interval
        table <- 
          table[-c(base::which(table$timestamp > 
                                 (table$timestamp[1] + seconds))),]
      }
    }
    table$subjectid <- ""
    table$subjectid[1] <- id
    table$subjectid[2] <- recordstart
    table$subjectid[3] <- recordstop
    table$subjectid <- as.character(table$subjectid)
    ##create column labeled "week", and fill with week "one" if < 7 days from recordstart, or "two" if > 7 days
    table$week <- ifelse(table$timestamp <= (table$timestamp[1]+lubridate::days(7)), "1", "2")
    table$week <- as.character(table$week)
    #extract the studyid based on the name of the unformatted file
    table$studyid <- sub("[ab].*","",basename(files[f]))
    table$subjectid[1] <- sub("[ab].*","",basename(files[f])) ## remove identifier from the file, and label with the studyid
    #extract the cgm period (a or b), based on the name of the unformatted file
    table$cgmperiod <- sub("DIG[0-9][0-9]*","",basename(files[f]))
    table$cgmperiod <- sub("[^ab].*","", table$cgmperiod)    
    table <-table[,c("subjectid","timestamp","sglucose", "cgmperiod", "week", "studyid")]
    
    
    ## if week = TRUE, the CGM periods is divided into 7 days blocks for analysis, otherwise is kept as fortnights.
    if(week == TRUE) {
    table1 <- dplyr::filter(table, week == "1") ## create separate tables for each week of CGM 
    table2 <- dplyr::filter(table, week == "2")
    
    table1$week <-  ifelse(table1$cgmperiod == "a", "1", "3") ## label weeks chronologically 1,2,3,4
    table2$week <-  ifelse(table2$cgmperiod == "b", "4", "2")
    
    filename1 <- 
      base::paste(output,"/",tools::file_path_sans_ext(
        basename(files[f])),"1.csv",sep = "")
    filename2 <- 
      base::paste(output,"/",tools::file_path_sans_ext(
        basename(files[f])),"2.csv",sep = "")
    utils::write.csv(as.data.frame(table1),file = filename1,row.names = FALSE)
    utils::write.csv(as.data.frame(table2),file = filename2,row.names = FALSE)
    } else {
      filename1 <- 
        base::paste(output,"/",tools::file_path_sans_ext(
          basename(files[f])),"1.csv",sep = "")
      utils::write.csv(as.data.frame(table),file = filename1,row.names = FALSE)
     
    }
   }}

## Run "tidycgm" function
tidycgm("Original", "Tidied", week = FALSE)

