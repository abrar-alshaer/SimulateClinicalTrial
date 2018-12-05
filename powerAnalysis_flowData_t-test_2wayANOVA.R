setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/Raz_Grant Nov2018/")
rm(list=ls())
library(pwr)

#determine effect sizes from previous study
data_t <- read.csv("EffectSizeData_lean-obese_flow_cell_studies.csv", header = TRUE) 
data_t_NonOb <- data_t[c(1:6),1]

#calculating effect size for 5S-HETE/14-HDHA study using Cohen's d method
#d = mean(pre) - mean(post) / pooled SD 
#instead of pre and post it can just be mean 1 and mean 2
d_data_t <- (mean(data_t_NonOb) - mean(data_t$Obese)) / sqrt(((sd(data_t_NonOb)^2 + sd(data_t$Obese)^2)/2))

pwr.t.test(n = 9, d = 1.45, sig.level = 0.05, type = c("two.sample"))
#You can specify alternative="two.sided", "less", or "greater" to indicate a two-tailed, or one-tailed test. A two tailed test is the default

#power analysis for 2-way ANOVA
#A & B are # of groups in factors A and B. We have 2 groups (DHA & control) and male and female for each of the 2 groups. 
#So we have 2 factors (M & F) with 2 groups each.
#beta = 1-power 

#cohen's d of 1.45 equals a Cohen's f of 0.7250 & eta-squared of 0.3445 -> converted the values from excel sheet "Converting effect sizes..."
#excel sheet web address: http://www.stat-help.com/spreadsheets/Converting%20effect%20sizes%202012-06-19.xls (http://www.stat-help.com/spreadsheets.html)
ss.2way(a=2,b=2,alpha=0.05,beta=0.1,f.A = 0.7250, f.B = 0.6, B=100)
#we are making the assumption here that both sex difference effect sizes & lean/obese DHA/control effect sizes are large

#Link with multiple power analyses: https://www.statmethods.net/stats/power.html 