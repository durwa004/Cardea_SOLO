###CS - import as .csv file
###TV - need to delete RR intervals by 1000 so that they are comparable

###First 24 hours - align CS and TV start points at the top of the spreadsheet
###Delete remaining CS values (below end of TV)

###Last 24 hours - align CS and TV end points at the end of the spreadsheet
###Delete remaining CS values (above start of TV)

#Read in RR intervals
RRs <- read.table("/Users/durwa004/Downloads/Nora_CS_TV_RR2.txt", header = T)
ts1 <- ts(RRs$CS_HR, start=1, frequency = 1)
ts2 <- ts(RRs$TV_HR, start=1, frequency = 1)
ccf_result <- ccf(ts1, ts2, lag.max = 10000) #Set for 30 minutes (also try for 60 minutes to be sure)
peak_lag <- ccf_result$lag[which.min(ccf_result$acf)]
peak_lag
min(ccf_result$acf)

