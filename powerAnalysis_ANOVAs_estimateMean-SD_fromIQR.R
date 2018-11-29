setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/Simulation Power Ananlysis/")
rm(list=ls())
library(pwr)
library(pwr2)
library(easypower)

#estimating mean & standard deviation from the median & IQR. 
#This method was adapted from the paper titled by: Wan, Xiang, Wenqian Wang, Jiming Liu, and Tiejun Tong. 2014. 
#"Estimating the Sample Mean and Standard Deviation from the Sample Size, Median, Range And/or Interquartile Range." 
#BMC Medical Research Methodology
#https://stats.stackexchange.com/questions/256456/how-to-calculate-mean-and-sd-from-quartiles 

#Estimating mean & SD for 30 days post-vaccination
q1 <- 70
median <- 170
q3 <- 200
n <- 40

mean_post <- (q1+median+q3)/3 #average of Q1,Q2,Q3
#you can only use the formula below if you have the sample size, if you don't have sample size use sd = q3-q1/1.35
sd_post <- (q3 - q1) / (2 * (qnorm((0.75 * n - 0.125) / (n + 0.25)))) 
#qnorm is used to estimate the the upper zth percentile of the standard normal distribution

#Estimating mean & SD for pre-vaccination
q1_pre <- 20
med_pre <- 25
q3_pre <- 60
n_pre <- 40
mean_pre <- (q1_pre+med_pre+q3_pre)/3
sd_pre <- (q3_pre - q1_pre) / (2 * (qnorm((0.75 * n_pre - 0.125) / (n_pre + 0.25)))) 

#calculating effect size for detecting HAI titers 30 days post vaccination
power <- (mean_pre - mean_post) / sqrt(((sd_pre^2 + sd_post^2)/2))

#power analysis for 2-way ANOVA
#A & B are # of groups in factors A and B. We have 2 groups (DHA & control) and male and female for each of the 2 groups. 
#So we have 2 factors (M & F) with 2 groups each. 
ss.2way(a=3,b=3,alpha=0.05,beta=0.2,f.A = 0.25, f.B = 0.25, B=100)
#Cohen's d 2-way ANOVA effect sizes:
#0.1 small effect size 
#0.25 medium effect size
#0.4 large effect size

#3-way ANOVA  
#all sample sizes are always rounded UP!
main.eff1 <- list(name = "Sex", levels = 2, eta.sq = "med")
main.eff2 <- list(name = "Time", levels = 3, eta.sq = "med")
main.eff3 <- list(name = "DHA", levels = 2, eta.sq = "med")
main.eff4 <- list(name = "BMI", levels = 3, eta.sq = "med")
n.multiway(iv1 = main.eff1, iv2 = main.eff2, iv3 = main.eff3, interaction.eta2 = "med", result = "all") #med is = 0.06 which is assuming a moderate effect size (med) for the interaction terms 
#you can specify the power you want with power = 0.90 or any other number. 
#4-way ANOVA
n.multiway(iv1 = main.eff1, iv2 = main.eff2, iv3 = main.eff3, iv4 = main.eff4, interaction.eta2 = 0.04, result = "all") #0.04 is estimating a slightly below moderate effect size for the interaction terms

#Link that contains the 3 way ANOVA code: https://cran.r-project.org/web/packages/easypower/vignettes/User_Input.html 

