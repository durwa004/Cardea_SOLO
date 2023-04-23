setwd("/Users/durwa004/Library/CloudStorage/Box-Box/Research/Projects/Cardiac/cardea_SOLO/Analysis/Animals_analysis/")
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
library(nlme)

#Get correlation/association between CS calculated HR and CS report HR
CS_HR <- read.table("CS_HR.txt", header=T)
CS_HR$dif <- CS_HR$CS_calc_HR - CS_HR$CS_rep_HR
cor.test(CS_HR$CS_calc_HR, CS_HR$CS_rep_HR)

#Import file and transform data
echos <- read.table("TD_combined_time.txt", header=T)
echos$Anon <- as.factor(echos$Anon)
echos$Time <- as.factor(echos$Time)

####BA plots of agreement between 2 methods
#Agreement between devices
#Can still do the BA analysis between the 2 devices - just can't do it between the raters
#Details: https://cran.r-project.org/web/packages/BlandAltmanLeh/vignettes/Intro.html

#HR
echos$HR_means <- (echos$TV_HR + echos$CS_HR)/2
echos$HR_diffs <- (echos$TV_HR - echos$CS_HR)
mean(echos$HR_diffs)
range(echos$HR_diffs)
mean(echos$HR_diffs) - (1.96*sd(echos$HR_diffs))
mean(echos$HR_diffs) + (1.96*sd(echos$HR_diffs))                      
echos$test <- echos$TV_HR
echos$CS <- echos$CS_HR
echos$max <- echos$test + (0.1*echos$test)
echos$min <- echos$test - (0.1*echos$test)
sum((echos$CS > echos$min) & (echos$CS < echos$max))/length(echos$test)

mean(echos$TV_HR[echos$Time == "A"])
mean(echos$TV_HR[echos$Time == "B"])
mean(echos$CS_HR[echos$Time == "A"])
mean(echos$CS_HR[echos$Time == "B"])

#ggplot
HR1.plot <- ggplot(echos,aes(x=HR_means, y = HR_diffs)) + 
  geom_point(alpha = 0.5, shape = echos$Time) + #circles = A, triangles = B
  geom_hline(yintercept = mean(echos$HR_diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(echos$HR_diffs) - (1.96*sd(echos$HR_diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(echos$HR_diffs) + (1.96*sd(echos$HR_diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO Heart rate (bpm") + 
  xlab("Average heart rate (bpm)") + 
  theme(panel.background = element_blank(), 
        plot.background = element_rect(colour = "NA"), 
        legend.position = "none", panel.grid = element_blank(), 
        panel.border = element_blank(),axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))
#HR
summary(lm(TV_HR ~ CS_HR + Time, data = echos))
HR1_m <- (lm(TV_HR ~ CS_HR + Time, data = echos))
confint(HR1_m)
echos$HR1_res <- HR1_m$residuals
ggplot(data = echos, aes(y = HR1_res, x = CS_HR)) + 
  geom_point(col = 'blue') + geom_abline(slope = 0)
bptest(HR1_m)


#Get EMMEANs
HR1_emm <- emmeans(HR1_m, specs = "Time", weights = "equal")
HR1_emm
x1 <- plot(HR1_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of heart rate (bpm)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#RR
echos$RR_means <- (echos$TV_RR + echos$CS_RR)/2
echos$RR_diffs <- (echos$TV_RR - echos$CS_RR)
mean(echos$RR_diffs) - (1.96*sd(echos$RR_diffs))
mean(echos$RR_diffs) + (1.96*sd(echos$RR_diffs))
mean(echos$RR_diffs)
range(echos$RR_diffs)
echos$test <- echos$TV_RR
echos$CS <- echos$CS_RR
echos$max <- echos$test + (0.1*echos$test)
echos$min <- echos$test - (0.1*echos$test)
sum((echos$CS > echos$min) & (echos$CS < echos$max), na.rm=T)/length(echos$test)

mean(echos$TV_RR[echos$Time == "A"])
mean(echos$TV_RR[echos$Time == "B"])
mean(echos$CS_RR[echos$Time == "A"])
mean(echos$CS_RR[echos$Time == "B"])

#ggplot
RR1.plot <- ggplot(echos,aes(x=RR_means, y = RR_diffs)) + 
  geom_point(alpha = 0.5, shape = echos$Time) +
  geom_hline(yintercept = mean(echos$RR_diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(echos$RR_diffs) - (1.96*sd(echos$RR_diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(echos$RR_diffs) + (1.96*sd(echos$RR_diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO RR interval (s)") + 
  xlab("Average RR interval (s)") + 
  theme(panel.background = element_blank(), 
        plot.background = element_rect(colour = "NA"), 
        legend.position = "none", panel.grid = element_blank(), 
        panel.border = element_blank(),axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))
#RR
summary(lm(TV_RR ~ CS_RR + Time, data = echos))
RR1_m <- (lm(TV_RR ~ CS_RR + Time, data = echos))

echos$RR1_res <- RR1_m$residuals
ggplot(data = echos, aes(y = RR1_res, x = CS_RR)) + 
  geom_point(col = 'blue') + geom_abline(slope = 0)
bptest(RR1_m)

#Get EMMEANs
RR1_emm <- emmeans(RR1_m, specs = "Time", weights = "equal")
x2 <- plot(RR1_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of RR interval (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#P_dur
echos$P_dur_means <- (echos$TV_P_dur + echos$CS_P_dur)/2
echos$P_dur_diffs <- (echos$TV_P_dur - echos$CS_P_dur)
mean(echos$P_dur_diffs, na.rm=T)
range(echos$P_dur_diffs)
mean(echos$P_dur_diffs) - (1.96*sd(echos$P_dur_diffs))
mean(echos$P_dur_diffs) + (1.96*sd(echos$P_dur_diffs))
echos$test <- echos$TV_P_dur
echos$CS <- echos$CS_P_dur
echos$max <- echos$test + (0.1*echos$test)
echos$min <- echos$test - (0.1*echos$test)
sum((echos$CS > echos$min) & (echos$CS < echos$max), na.rm=T)/length(echos$test)

mean(echos$TV_P_dur[echos$Time == "A"])
mean(echos$TV_P_dur[echos$Time == "B"])
mean(echos$CS_P_dur[echos$Time == "A"])
mean(echos$CS_P_dur[echos$Time == "B"])

#ggplot
P_dur1.plot <- ggplot(echos,aes(x=P_dur_means, y = P_dur_diffs)) + 
  geom_point(alpha = 0.5, shape = echos$Time) +
  geom_hline(yintercept = mean(echos$P_dur_diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(echos$P_dur_diffs) - (1.96*sd(echos$P_dur_diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(echos$P_dur_diffs) + (1.96*sd(echos$P_dur_diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO P duration (s)") + 
  xlab("Average P duration (s)") + 
  theme(panel.background = element_blank(), 
        plot.background = element_rect(colour = "NA"), 
        legend.position = "none", panel.grid = element_blank(), 
        panel.border = element_blank(),axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#P_dur
summary(lm(TV_P_dur ~ CS_P_dur + Time, data = echos))
P_dur1_m <- (lm(TV_P_dur ~ CS_P_dur + Time, data = echos))

echos$P_dur1_res <- P_dur1_m$residuals
ggplot(data = echos, aes(y = P_dur1_res, x = CS_P_dur)) + 
  geom_point(col = 'blue') + geom_abline(slope = 0)
bptest(P_dur1_m)

#Get EMMEANs
P_dur1_emm <- emmeans(P_dur1_m, specs = "Time", weights = "equal")
x3 <- plot(P_dur1_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of P duration (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#PRI
echos$PRI_means <- (echos$TV_PRI + echos$CS_PRI)/2
echos$PRI_diffs <- (echos$TV_PRI - echos$CS_PRI)
mean(echos$PRI_diffs)
range(echos$PRI_diffs)
mean(echos$PRI_diffs) - (1.96*sd(echos$PRI_diffs))
mean(echos$PRI_diffs) + (1.96*sd(echos$PRI_diffs))
echos$test <- echos$TV_PRI
echos$CS <- echos$CS_PRI
echos$max <- echos$test + (0.1*echos$test)
echos$min <- echos$test - (0.1*echos$test)
sum((echos$CS > echos$min) & (echos$CS < echos$max), na.rm=T)/length(echos$test)

mean(echos$TV_PRI[echos$Time == "A"])
mean(echos$TV_PRI[echos$Time == "B"])
mean(echos$CS_PRI[echos$Time == "A"])
mean(echos$CS_PRI[echos$Time == "B"])

#ggplot
PRI1.plot <- ggplot(echos,aes(x=PRI_means, y = PRI_diffs)) + 
  geom_point(alpha = 0.5, shape = echos$Time) +
  geom_hline(yintercept = mean(echos$PRI_diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(echos$PRI_diffs) - (1.96*sd(echos$PRI_diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(echos$PRI_diffs) + (1.96*sd(echos$PRI_diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO PR interval (s)") + 
  xlab("Average PR interval (s)") + 
  theme(panel.background = element_blank(), 
        plot.background = element_rect(colour = "NA"), 
        legend.position = "none", panel.grid = element_blank(), 
        panel.border = element_blank(),axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))
#PRI
summary(lm(TV_PRI ~ CS_PRI + Time, data = echos))
PRI1_m <- (lm(TV_PRI ~ CS_PRI + Time, data = echos))

echos$PRI1_res <- PRI1_m$residuals
ggplot(data = echos, aes(y = PRI1_res, x = CS_PRI)) + 
  geom_point(col = 'blue') + geom_abline(slope = 0)
bptest(PRI1_m)

#Get EMMEANs
PRI1_emm <- emmeans(PRI1_m, specs = "Time", weights = "equal")
x4 <- plot(PRI1_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of PR interval (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#QTI
echos$QTI_means <- (echos$TV_QTI + echos$CS_QTI)/2
echos$QTI_diffs <- (echos$TV_QTI - echos$CS_QTI)
mean(echos$QTI_diffs, na.rm=T)
range(echos$QTI_diffs)
mean(echos$QTI_diffs) - (1.96*sd(echos$QTI_diffs))
mean(echos$QTI_diffs) + (1.96*sd(echos$QTI_diffs))
echos$test <- echos$TV_QTI
echos$CS <- echos$CS_QTI
echos$max <- echos$test + (0.1*echos$test)
echos$min <- echos$test - (0.1*echos$test)
sum((echos$CS > echos$min) & (echos$CS < echos$max), na.rm=T)/length(echos$test)

mean(echos$TV_QTI[echos$Time == "A"])
mean(echos$TV_QTI[echos$Time == "B"])
mean(echos$CS_QTI[echos$Time == "A"])
mean(echos$CS_QTI[echos$Time == "B"])

#ggplot
QTI1.plot <- ggplot(echos,aes(x=QTI_means, y = QTI_diffs)) + 
  geom_point(alpha = 0.5, shape = echos$Time) +
  geom_hline(yintercept = mean(echos$QTI_diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(echos$QTI_diffs) - (1.96*sd(echos$QTI_diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(echos$QTI_diffs) + (1.96*sd(echos$QTI_diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO QT interval (s)") + 
  xlab("Average QT interval (s)") + 
  theme(panel.background = element_blank(), 
        plot.background = element_rect(colour = "NA"), 
        legend.position = "none", panel.grid = element_blank(), 
        panel.border = element_blank(),axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))
#QTI
summary(lm(TV_QTI ~ CS_QTI + Time, data = echos))
QTI1_m <- (lm(TV_QTI ~ CS_QTI + Time, data = echos))

echos$QTI1_res <- QTI1_m$residuals
ggplot(data = echos, aes(y = QTI1_res, x = CS_QTI)) + 
  geom_point(col = 'blue') + geom_abline(slope = 0)
bptest(QTI1_m)
QTI1_m2 <- nlme::gls(TV_QTI ~ CS_QTI + Time, 
                     weights = varIdent( ~1 | CS_QTI), data = echos)
#Get EMMEANs
QTI1_emm <- emmeans(QTI1_m2, specs = "Time", weights = "equal")
x5 <- plot(PRI1_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of QT interval (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#PRS
echos$PRS_means <- (echos$TV_PRS + echos$CS_PRS)/2
echos$PRS_diffs <- (echos$TV_PRS - echos$CS_PRS)
mean(echos$PRS_diffs, na.rm=T)
range(echos$PRS_diffs)
mean(echos$PRS_diffs) - (1.96*sd(echos$PRS_diffs))
mean(echos$PRS_diffs) + (1.96*sd(echos$PRS_diffs))
echos$test <- echos$TV_PRS
echos$CS <- echos$CS_PRS
echos$max <- echos$test + (0.1*echos$test)
echos$min <- echos$test - (0.1*echos$test)
sum((echos$CS > echos$min) & (echos$CS < echos$max), na.rm=T)/length(echos$test)

mean(echos$TV_PRS[echos$Time == "A"])
mean(echos$TV_PRS[echos$Time == "B"])
mean(echos$CS_PRS[echos$Time == "A"])
mean(echos$CS_PRS[echos$Time == "B"])

#ggplot
PRS1.plot <- ggplot(echos,aes(x=PRS_means, y = PRS_diffs)) + 
  geom_point(alpha = 0.5, shape = echos$Time) +
  geom_hline(yintercept = mean(echos$PRS_diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(echos$PRS_diffs) - (1.96*sd(echos$PRS_diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(echos$PRS_diffs) + (1.96*sd(echos$PRS_diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO PR segment (s)") + 
  xlab("Average PR segment (s)") + 
  theme(panel.background = element_blank(), 
        plot.background = element_rect(colour = "NA"), 
        legend.position = "none", panel.grid = element_blank(), 
        panel.border = element_blank(),axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))
#PRS
summary(lm(TV_PRS ~ CS_PRS + Time, data = echos))
PRS1_m <- (lm(TV_PRS ~ CS_PRS + Time, data = echos))

echos$PRS1_res <- PRS1_m$residuals
ggplot(data = echos, aes(y = PRS1_res, x = CS_PRS)) + 
  geom_point(col = 'blue') + geom_abline(slope = 0)
bptest(PRS1_m)
PRS1_m2 <- nlme::gls(TV_PRS ~ CS_PRS + Time, 
               weights = varIdent( ~1 | CS_PRS), data = echos)

#Get EMMEANs
PRS1_emm <- emmeans(PRS1_m2, specs = "Time", weights = "equal")
x6 <- plot(PRS1_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of PR segment (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#STS
echos$STS_means <- (echos$TV_STS + echos$CS_STS)/2
echos$STS_diffs <- (echos$TV_STS - echos$CS_STS)
mean(echos$STS_diffs, na.rm=T)
range(echos$STS_diffs)
mean(echos$STS_diffs) - (1.96*sd(echos$STS_diffs))
mean(echos$STS_diffs) + (1.96*sd(echos$STS_diffs))
echos$test <- echos$TV_STS
echos$CS <- echos$CS_STS
echos$max <- echos$test + (0.1*echos$test)
echos$min <- echos$test - (0.1*echos$test)
sum((echos$CS > echos$min) & (echos$CS < echos$max), na.rm=T)/length(echos$test)

mean(echos$TV_STS[echos$Time == "A"])
mean(echos$TV_STS[echos$Time == "B"])
mean(echos$CS_STS[echos$Time == "A"])
mean(echos$CS_STS[echos$Time == "B"])

#ggplot
STS1.plot <- ggplot(echos,aes(x=STS_means, y = STS_diffs)) + 
  geom_point(alpha = 0.5, shape = echos$Time) +
  geom_hline(yintercept = mean(echos$STS_diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(echos$STS_diffs) - (1.96*sd(echos$STS_diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(echos$STS_diffs) + (1.96*sd(echos$STS_diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO ST segment (s)") + 
  xlab("Average ST segment (s)") + 
  theme(panel.background = element_blank(), 
        plot.background = element_rect(colour = "NA"), 
        legend.position = "none", panel.grid = element_blank(), 
        panel.border = element_blank(),axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))
#STS
summary(lm(TV_STS ~ CS_STS + Time, data = echos))
STS1_m <- (lm(TV_STS ~ CS_STS + Time, data = echos))

echos$STS1_res <- STS1_m$residuals
ggplot(data = echos, aes(y = STS1_res, x = CS_STS)) + 
  geom_point(col = 'blue') + geom_abline(slope = 0)
bptest(STS1_m)

#Get EMMEANs
STS1_emm <- emmeans(STS1_m, specs = "Time", weights = "equal")
x7 <- plot(STS1_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of ST segment (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

#QRS_dur
echos$QRS_dur_means <- (echos$TV_QRS_dur + echos$CS_QRS_dur)/2
echos$QRS_dur_diffs <- (echos$TV_QRS_dur - echos$CS_QRS_dur)
mean(echos$QRS_dur_diffs, na.rm=T)
range(echos$QRS_dur_diffs)
mean(echos$QRS_dur_diffs) - (1.96*sd(echos$QRS_dur_diffs))
mean(echos$QRS_dur_diffs) + (1.96*sd(echos$QRS_dur_diffs))
echos$test <- echos$TV_QRS_dur
echos$CS <- echos$CS_QRS_dur
echos$max <- echos$test + (0.1*echos$test)
echos$min <- echos$test - (0.1*echos$test)
sum((echos$CS > echos$min) & (echos$CS < echos$max), na.rm=T)/length(echos$test)

mean(echos$TV_QRS_dur[echos$Time == "A"])
mean(echos$TV_QRS_dur[echos$Time == "B"])
mean(echos$CS_QRS_dur[echos$Time == "A"])
mean(echos$CS_QRS_dur[echos$Time == "B"])

#ggplot
QRS_dur1.plot <- ggplot(echos,aes(x=QRS_dur_means, y = QRS_dur_diffs)) + 
  geom_point(alpha = 0.5, shape = echos$Time) +
  geom_hline(yintercept = mean(echos$QRS_dur_diffs), colour = "blue", size = 0.5) + 
  geom_hline(yintercept = mean(echos$QRS_dur_diffs) - (1.96*sd(echos$QRS_dur_diffs)), colour = "red", size = 0.5, linetype = 3) +
  geom_hline(yintercept = mean(echos$QRS_dur_diffs) + (1.96*sd(echos$QRS_dur_diffs)), colour = "red", size = 0.5, linetype = 3) + 
  ylab("Difference between Televet and\nCardea SOLO QRS duration (s)") + 
  xlab("Average QRS duration (s)") + 
  theme(panel.background = element_blank(), 
        plot.background = element_rect(colour = "NA"), 
        legend.position = "none", panel.grid = element_blank(), 
        panel.border = element_blank(),axis.line.x = element_line(),
        axis.line.y = element_line(), axis.text.x = element_text(), 
        axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))
#QRS_dur
summary(lm(TV_QRS_dur ~ CS_QRS_dur + Time, data = echos))
QRS_dur1_m <- (lm(TV_QRS_dur ~ CS_QRS_dur + Time, data = echos))

echos$QRS_dur1_res <- QRS_dur1_m$residuals
ggplot(data = echos, aes(y = QRS_dur1_res, x = CS_QRS_dur)) + 
  geom_point(col = 'blue') + geom_abline(slope = 0)
bptest(QRS_dur1_m)


#Get EMMEANs
QRS_dur1_emm <- emmeans(QRS_dur1_m, specs = "Time", weights = "equal")
x8 <- plot(QRS_dur1_emm) + geom_boxplot() + theme_bw() + xlab("EMMEAN of QRS duraton (s)") + 
  ylab("Time period") +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.line.x = element_line(), axis.line.y = element_line(), 
        axis.text.x = element_text(), axis.text = element_text(size=10), 
        axis.title = element_text(size=12,face="bold"))

BA <- plot_grid(HR1.plot, P_dur1.plot, QRS_dur1.plot, PRI1.plot, RR1.plot, QTI1.plot, PRS1.plot, STS1.plot, 
                  labels = c("A","B","C","D","E", "F", "G", "H"),
                  ncol = 2)
save_plot("BA.pdf", BA,ncol=2,nrow=4)

EM <- plot_grid(x1, x3, x8, x4, x2,x5, x6, x7, 
                labels = c("A","B","C","D","E", "F", "G", "H"),
                ncol = 2)
save_plot("EM.pdf", EM,ncol=2,nrow=4)

##Agreement: We can use the ICC(1,1) model to compare two measurements of the same variable where one is considered to be the “gold standard”.
icc(echos[,c("TV1_HR","CS1_HR")], model = "twoway", type = "consistency")
icc(echos[,c("TV1_RR","CS1_RR")], model = "twoway", type = "consistency")
icc(echos[,c("TV1_P_dur","CS1_P_dur")], model = "twoway", type = "consistency")
icc(echos[,c("TV1_PRI","CS1_PRI")], model = "twoway", type = "consistency")
icc(echos[,c("TV1_QTI","CS1_QTI")], model = "twoway", type = "consistency")
icc(echos[,c("TV1_PRS","CS1_PRS")], model = "twoway", type = "consistency")
icc(echos[,c("TV1_STS","CS1_STS")], model = "twoway", type = "consistency")
icc(echos[,c("TV1_QRS_dur","CS1_QRS_dur")], model = "twoway", type = "consistency")

icc(echos[,c("TV2_HR","CS2_HR")], model = "twoway", type = "consistency")
icc(echos[,c("TV2_RR","CS2_RR")], model = "twoway", type = "consistency")
icc(echos[,c("TV2_P_dur","CS2_P_dur")], model = "twoway", type = "consistency")
icc(echos[,c("TV2_PRI","CS2_PRI")], model = "twoway", type = "consistency")
icc(echos[,c("TV2_QTI","CS2_QTI")], model = "twoway", type = "consistency")
icc(echos[,c("TV2_PRS","CS2_PRS")], model = "twoway", type = "consistency")
icc(echos[,c("TV2_STS","CS2_STS")], model = "twoway", type = "consistency")
icc(echos[,c("TV2_QRS_dur","CS2_QRS_dur")], model = "twoway", type = "consistency")

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

