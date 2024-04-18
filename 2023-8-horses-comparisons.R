setwd("/Users/durwa004/Library/CloudStorage/Box-Box/Research/Projects/Cardiac/cardea_SOLO/Analysis/8_horses_comparison/")
##Correlation = correlation coefficients
#Check agreement between AC and TV Cohen's kappa and using BA plots
library(BlandAltmanLeh)
library(ggplot2)
library(cowplot)
library(irr)
library(dplyr)
library(forcats)
library(devtools)
library(ggpubr)
library(emmeans)
library(scales)
library(dvmisc)
library(tidyr)
library(lmtest)

#Import file and transform data
echos <- read.table("8_horses_combined.txt", header=T)
#echos$TV_arrhythmia <- as.factor(echos$TV_arrhythmia)
#echos$TV_PAC <- as.factor(echos$TV_PAC)
#echos$TV_PVC <- as.factor(echos$TV_PVC)
#echos$CS_arrhythmia <- as.factor(echos$CS_arrhythmia)
#echos$CS_PAC <- as.factor(echos$CS_PAC)
#echos$CS_PVC <- as.factor(echos$CS_PVC)
echos$TV_HR <- rowMeans(echos[,c(2:4)], na.rm=T)
echos$TV_RR <- (rowMeans(echos[,c(5:7)], na.rm=T))/1000
echos$TV_P_dur <- rowMeans(echos[,c(8:10)], na.rm=T)
echos$TV_PRI <- rowMeans(echos[,c(11:13)], na.rm=T)
echos$TV_QTI <- rowMeans(echos[,c(14:16)], na.rm=T)
echos$TV_PRS <- rowMeans(echos[,c(17:19)], na.rm=T)
echos$TV_STS <- rowMeans(echos[,c(20:22)], na.rm=T)
echos$TV_QRS_dur <- rowMeans(echos[,c(23:25)], na.rm=T)

echos$CS_HR <- rowMeans(echos[,c(26:28)], na.rm=T)
echos$CS_RR <- (rowMeans(echos[,c(29:31)], na.rm=T))/1000
echos$CS_P_dur <- rowMeans(echos[,c(32:34)], na.rm=T)
echos$CS_PRI <- rowMeans(echos[,c(35:37)], na.rm=T)
echos$CS_QTI <- rowMeans(echos[,c(38:40)], na.rm=T)
echos$CS_PRS <- rowMeans(echos[,c(41:43)], na.rm=T)
echos$CS_STS <- rowMeans(echos[,c(44:46)], na.rm=T)
echos$CS_QRS_dur <- rowMeans(echos[,c(47:49)], na.rm=T)

##Correlations
cor(echos$TV_RR, echos$CS_RR, method = "pearson", use = "complete.obs")
cor(echos$TV_P_dur, echos$CS_P_dur, method = "pearson", use = "complete.obs")
cor(echos$TV_PRI, echos$CS_PRI, method = "pearson", use = "complete.obs")
cor(echos$TV_QTI, echos$CS_QTI, method = "pearson", use = "complete.obs")
cor(echos$TV_PRS, echos$CS_PRS, method = "pearson", use = "complete.obs")
cor(echos$TV_STS, echos$CS_STS, method = "pearson", use = "complete.obs")
cor(echos$TV_QRS_dur, echos$CS_QRS_dur, method = "pearson", use = "complete.obs")


#% of Cardea SOLO values within 10% of TV value
echos$test <- echos$TV_RR
echos$CS <- echos$CS_RR
echos$max <- echos$test + (0.1*echos$test)
echos$min <- echos$test - (0.1*echos$test)

(sum((echos$CS > echos$min) & (echos$CS < echos$max), na.rm=T)/length(echos$test))*100




####BA plots of agreement between 2 methods
#Agreement between devices
#Can still do the BA analysis between the 2 devices - just can't do it between the raters
#Details: https://cran.r-project.org/web/packages/BlandAltmanLeh/vignettes/Intro.html

#HR
#Use .stats file to get median paired difference
p.stats <- bland.altman.stats(echos$TV_HR,echos$CS_HR)
mean(p.stats$diffs)
df.p <- data.frame(p.stats$means,p.stats$diffs)
#ggplot
HR1.plot <- ggplot(df.p,aes(x=p.stats.means, y = p.stats.diffs)) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs), colour = "blue", linewidth = 0.5) + 
  geom_hline(yintercept = mean(df.p$p.stats.diffs) - (1.96*sd(df.p$p.stats.diffs)), colour = "red", linewidth = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs) + (1.96*sd(df.p$p.stats.diffs)), colour = "red", linewidth = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO Heart rate") + 
  xlab("Average heart rate (bpm)") + theme(panel.background = element_blank(), 
                                         plot.background = element_rect(colour = "NA"), legend.position = "none", 
                                         panel.grid = element_blank(), panel.border = element_blank(),axis.line.x = element_line(),
                                         axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
summary(df.p$p.stats.diffs)
#RR
#Use .stats file to get median paired difference
p.stats <- bland.altman.stats(echos$TV_RR,echos$CS_RR)
mean(p.stats$diffs) #0.20
#ggplot
df.p <- data.frame(p.stats$means,p.stats$diffs)
RR1.plot <- ggplot(df.p,aes(x=p.stats.means, y = p.stats.diffs)) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(df.p$p.stats.diffs) - (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs) + (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO RR interval") + 
  xlab("Average RR interval (ms)") + theme(panel.background = element_blank(), 
                                     plot.background = element_rect(colour = "NA"), legend.position = "none", 
                                     panel.grid = element_blank(), panel.border = element_blank(),axis.line.x = element_line(),
                                     axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

summary(df.p$p.stats.diffs)
#P_dur
#Use .stats file to get median paired difference
p.stats <- bland.altman.stats(echos$TV_P_dur,echos$CS_P_dur)
mean(p.stats$diffs) #0.20
#ggplot
df.p <- data.frame(p.stats$means,p.stats$diffs)
P_dur1.plot <- ggplot(df.p,aes(x=p.stats.means, y = p.stats.diffs)) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(df.p$p.stats.diffs) - (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs) + (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO P duration") + 
  xlab("Average P duration (s)") + theme(panel.background = element_blank(), 
                                          plot.background = element_rect(colour = "NA"), legend.position = "none", 
                                          panel.grid = element_blank(), panel.border = element_blank(),axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
summary(df.p$p.stats.diffs)
#PRI
#Use .stats file to get median paired difference
p.stats <- bland.altman.stats(echos$TV_PRI,echos$CS_PRI)
mean(p.stats$diffs) #0.20
#ggplot
df.p <- data.frame(p.stats$means,p.stats$diffs)
PRI1.plot <- ggplot(df.p,aes(x=p.stats.means, y = p.stats.diffs)) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(df.p$p.stats.diffs) - (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs) + (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO PR interval") + 
  xlab("Average PR interval (s)") + theme(panel.background = element_blank(), 
                                          plot.background = element_rect(colour = "NA"), legend.position = "none", 
                                          panel.grid = element_blank(), panel.border = element_blank(),axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
summary(df.p$p.stats.diffs)

#QTI
#Use .stats file to get median paired difference
p.stats <- bland.altman.stats(echos$TV_QTI,echos$CS_QTI)
mean(p.stats$diffs) #0.20
#ggplot
df.p <- data.frame(p.stats$means,p.stats$diffs)
QTI1.plot <- ggplot(df.p,aes(x=p.stats.means, y = p.stats.diffs)) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(df.p$p.stats.diffs) - (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs) + (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO QT interval") + 
  xlab("Average QT interval (s)") + theme(panel.background = element_blank(), 
                                          plot.background = element_rect(colour = "NA"), legend.position = "none", 
                                          panel.grid = element_blank(), panel.border = element_blank(),axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
summary(df.p$p.stats.diffs)

#PRS
#Use .stats file to get median paired difference
p.stats <- bland.altman.stats(echos$TV_PRS,echos$CS_PRS)
mean(p.stats$diffs) #0.20
#ggplot
df.p <- data.frame(p.stats$means,p.stats$diffs)
PRS1.plot <- ggplot(df.p,aes(x=p.stats.means, y = p.stats.diffs)) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(df.p$p.stats.diffs) - (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs) + (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO PR segment duration") + 
  xlab("Average PR segment duration (s)") + theme(panel.background = element_blank(), 
                                          plot.background = element_rect(colour = "NA"), legend.position = "none", 
                                          panel.grid = element_blank(), panel.border = element_blank(),axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
summary(df.p$p.stats.diffs)

#STS
#Use .stats file to get median paired difference
p.stats <- bland.altman.stats(echos$TV_STS,echos$CS_STS)
mean(p.stats$diffs)
#ggplot
df.p <- data.frame(p.stats$means,p.stats$diffs)
STS1.plot <- ggplot(df.p,aes(x=p.stats.means, y = p.stats.diffs)) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(df.p$p.stats.diffs) - (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs) + (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO ST segment duration") + 
  xlab("Average ST segment duration (s)") + theme(panel.background = element_blank(), 
                                          plot.background = element_rect(colour = "NA"), legend.position = "none", 
                                          panel.grid = element_blank(), panel.border = element_blank(),axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

summary(df.p$p.stats.diffs)

#QRS_dur
#Use .stats file to get median paired difference
p.stats <- bland.altman.stats(echos$TV_QRS_dur,echos$CS_QRS_dur)
mean(p.stats$diffs)
#ggplot
df.p <- data.frame(p.stats$means,p.stats$diffs)
QRS_dur1.plot <- ggplot(df.p,aes(x=p.stats.means, y = p.stats.diffs)) + 
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(df.p$p.stats.diffs) - (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(df.p$p.stats.diffs) + (1.96*sd(df.p$p.stats.diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO QRS duration") + 
  xlab("Average QRS duration (s)") + theme(panel.background = element_blank(), 
                                          plot.background = element_rect(colour = "NA"), legend.position = "none", 
                                          panel.grid = element_blank(), panel.border = element_blank(),axis.line.x = element_line(),
                                          axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))
summary(df.p$p.stats.diffs)

BA.1 <- plot_grid(HR1.plot, P_dur1.plot, QRS_dur1.plot, PRI1.plot, RR1.plot, QTI1.plot, PRS1.plot, STS1.plot, 
                    labels = c("A","B","C","D","E", "F", "G", "H"),
                    ncol = 2)
save_plot("BA_8_horse_comparison.pdf", BA.1,ncol=2,nrow=4)


##Agreement: We can use the ICC(1,1) model to compare two measurements of the same variable where one is considered to be the “gold standard”.
icc(echos[,c("TV_HR","CS_HR")], model = "twoway", type = "consistency")
icc(echos[,c("TV_RR","CS_RR")], model = "twoway", type = "consistency")
icc(echos[,c("TV_P_dur","CS_P_dur")], model = "twoway", type = "consistency")
icc(echos[,c("TV_PRI","CS_PRI")], model = "twoway", type = "consistency")
icc(echos[,c("TV_QTI","CS_QTI")], model = "twoway", type = "consistency")
icc(echos[,c("TV_PRS","CS_PRS")], model = "twoway", type = "consistency")
icc(echos[,c("TV_STS","CS_STS")], model = "twoway", type = "consistency")
icc(echos[,c("TV_QRS_dur","CS_QRS_dur")], model = "twoway", type = "consistency")

#Usability of cardea solo data
CS_data <- read.table("Hourly_readability.txt", header=T)
sum(rowSums(CS_data == "NSR", na.rm=T))/(sum(sum(rowSums(CS_data == "NSR", na.rm=T))) + 
                                           sum(rowSums(CS_data == "ST", na.rm=T)) +
                                           sum(rowSums(CS_data == "None", na.rm=T)) +
                                           sum(rowSums(CS_data == "UR", na.rm=T)))
sum(rowSums(CS_data == "ST", na.rm=T))/(sum(sum(rowSums(CS_data == "NSR", na.rm=T))) + 
                                           sum(rowSums(CS_data == "ST", na.rm=T)) +
                                           sum(rowSums(CS_data == "None", na.rm=T)) +
                                           sum(rowSums(CS_data == "UR", na.rm=T)))
sum(rowSums(CS_data == "UR", na.rm=T))/(sum(sum(rowSums(CS_data == "NSR", na.rm=T))) + 
                                           sum(rowSums(CS_data == "ST", na.rm=T)) +
                                           sum(rowSums(CS_data == "None", na.rm=T)) +
                                           sum(rowSums(CS_data == "UR", na.rm=T)))
sum(rowSums(CS_data == "None", na.rm=T))/(sum(sum(rowSums(CS_data == "NSR", na.rm=T))) + 
                                           sum(rowSums(CS_data == "ST", na.rm=T)) +
                                           sum(rowSums(CS_data == "None", na.rm=T)) +
                                           sum(rowSums(CS_data == "UR", na.rm=T)))

CS_data$NSR <- rowSums(CS_data == "NSR", na.rm=T)
CS_data$ST <- rowSums(CS_data == "ST", na.rm=T)
CS_data$UR <- rowSums(CS_data == "UR", na.rm=T)
CS_data$None <- rowSums(CS_data == "None", na.rm=T)

CS_data$readable <- CS_data$NSR + CS_data$ST
CS_data$unreadable <- CS_data$UR + CS_data$None


day1 <- c(sum(CS_data$readable[CS_data$Day == 1]), sum(CS_data$unreadable[CS_data$Day == 1]))
day2 <- c(sum(CS_data$readable[CS_data$Day == 2]), sum(CS_data$unreadable[CS_data$Day == 2]))
day3 <- c(sum(CS_data$readable[CS_data$Day == 3]), sum(CS_data$unreadable[CS_data$Day == 3]))
day4 <- c(sum(CS_data$readable[CS_data$Day == 4]), sum(CS_data$unreadable[CS_data$Day == 4]))
day5 <- c(sum(CS_data$readable[CS_data$Day == 5]), sum(CS_data$unreadable[CS_data$Day == 5]))
day6 <- c(sum(CS_data$readable[CS_data$Day == 6]), sum(CS_data$unreadable[CS_data$Day == 6]))
day7 <- c(sum(CS_data$readable[CS_data$Day == 7]), sum(CS_data$unreadable[CS_data$Day == 7]))

CS_brief <- as.data.frame(rbind(day1,day7))
chisq.test(CS_brief)

day <- (c(rep("1",2), rep("2",2), rep("3",2), rep("4",2), 
          rep("5",2), rep("6", 2), rep("7", 2)))
condition <- rep(c("readable", "unreadable"), 7)
results <- c(67,5, 71,1,70,2,62,10,55,17,39,33, 35, 24)

CS_plot <- as.data.frame(day, condition)
CS_plot$results <- as.numeric(results)

ECG.plot <- ggplot(CS_plot, aes(fill=condition, y=results, x=day)) + 
  theme_bw() + 
  geom_bar(position="dodge", stat="identity") + 
  ylab("Count") + 
  xlab("Day") + scale_y_continuous(limits = c(0,80)) + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), 
        axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"),
        legend.title = element_blank(), legend.position = "none")
save_plot("ECG_quality.tiff", ECG.plot, base_height = 6, base_width = 6, dpi = 300)




##For animals review - don't think we need

#RR
summary(lm(TV_RR ~ CS_RR + Time, data = echos))
RR1_m <- (lm(TV_RR ~ CS_RR + Time, data = echos))
bptest(RR1_m)

#Get EMMEANs
RR1_emm <- emmeans(RR1_m, specs = "Time", weights = "equal")
x2 <- plot(RR1_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of RR interval (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#P duration
summary(lm(TV_P_dur ~ CS_P_dur + Time, data = echos))
P_dur_m <- (lm(TV_P_dur ~ CS_P_dur + Time, data = echos))
bptest(P_dur_m)

#Get EMMEANs
P_dur_emm <- emmeans(P_dur_m, specs = "Time", weights = "equal")
x3 <- plot(P_dur_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of P duration (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#PRI
summary(lm(TV_PRI ~ CS_PRI + Time, data = echos))
PRI_m <- (lm(TV_PRI ~ CS_PRI + Time, data = echos))
bptest(PRI_m)

#Get EMMEANs
PRI_emm <- emmeans(PRI_m, specs = "Time", weights = "equal")
x4 <- plot(PRI_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of PR interval (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#QTI
summary(lm(TV_QTI ~ CS_QTI + Time, data = echos))
QTI_m <- (lm(TV_QTI ~ CS_QTI + Time, data = echos))
bptest(QTI_m)
#Get EMMEANs
QTI_emm <- emmeans(QTI_m, specs = "Time", weights = "equal")
x5 <- plot(QTI_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of QT interval (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#PRS
summary(lm(TV_PRS ~ CS_PRS + Time, data = echos))
PRS_m <- (lm(TV_PRS ~ CS_PRS + Time, data = echos))
bptest(PRS_m)

#Get EMMEANs
PRS_emm <- emmeans(PRS_m, specs = "Time", weights = "equal")
x6 <- plot(PRS_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of PR segment (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#STS
summary(lm(TV_STS ~ CS_STS + Time, data = echos))
STS_m <- (lm(TV_STS ~ CS_STS + Time, data = echos))
bptest(STS_m)

#Get EMMEANs
STS_emm <- emmeans(STS_m, specs = "Time", weights = "equal")
x7 <- plot(STS_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of ST segment (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

#QRS duration
summary(lm(TV_QRS_dur ~ CS_QRS_dur + Time, data = echos))
QRS_dur_m <- (lm(TV_QRS_dur ~ CS_QRS_dur + Time, data = echos))
bptest(QRS_dur_m)

#Get EMMEANs
QRS_dur_emm <- emmeans(QRS_dur_m, specs = "Time", weights = "equal")
x8 <- plot(QRS_dur_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of QRS duration (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), axis.text = element_text(size=10), axis.title = element_text(size=12,face="bold"))

