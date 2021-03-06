setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/Raz_Grant Nov2018/")
rm(list=ls())
library(pwr)
library(pwr2)
library(easypower)

#estimating mean & standard deviation from the median & IQR. 
#This method was adapted from the paper titled by: Wan, Xiang, Wenqian Wang, Jiming Liu, and Tiejun Tong. 2014. 
#“Estimating the Sample Mean and Standard Deviation from the Sample Size, Median, Range And/or Interquartile Range.” 
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

#power analysis for T-test
pwr.t.test(n = 228, d = 0.3, sig.level = 0.05, type = c("two.sample"))
#You can specify alternative="two.sided", "less", or "greater" to indicate a two-tailed, or one-tailed test. A two tailed test is the default

#power analysis for One-way ANOVA
pwr.anova.test(k = 2, n = NULL, f = 0.2, sig.level = 0.05, power = 0.8)
#where k is the number of groups and n is the common sample size in each group, f = effect size
#effect size options are the same as below.

#power analysis for 2-way ANOVA
#A & B are # of groups in factors A and B. We have 2 groups (DHA & control) and male and female for each of the 2 groups. 
#So we have 2 factors (M & F) with 2 groups each. 
ss.2way(a=2,b=2,alpha=0.05,beta=0.2,f.A = 0.2, f.B = 0.2, B=100)
#Cohen's d 2-way ANOVA effect sizes:
#0.1 small effect size 
#0.25 medium effect size
#0.4 large effect size

#3-way ANOVA  
#all sample sizes are always rounded UP!
main.eff1 <- list(name = "Sex", levels = 2, eta.sq = "med")
main.eff2 <- list(name = "SNPs", levels = 3, eta.sq = "med")
main.eff3 <- list(name = "DHA", levels = 2, eta.sq = "med")
main.eff4 <- list(name = "BMI", levels = 3, eta.sq = "med")
n.multiway(iv1 = main.eff1, iv2 = main.eff2, iv3 = main.eff3, interaction.eta2 = "med", result = "all", power = 0.80) #med is = 0.06 which is assuming a moderate effect size (med) for the interaction terms 
#you can specify the power you want with power = 0.90 or any other number. 
#4-way ANOVA
n.multiway(iv1 = main.eff1, iv2 = main.eff2, iv3 = main.eff3, iv4 = main.eff4, interaction.eta2 = 0.04, result = "all") #0.04 is estimating a slightly below moderate effect size for the interaction terms

#Link that contains the 3 way ANOVA code: https://cran.r-project.org/web/packages/easypower/vignettes/User_Input.html 

#calculating effect sizes for pearson & spearman correlation coefficiants
pwr.r.test(n = NULL, r = 0.7, sig.level = 0.05, power = 0.8) #n= NULL because we want an estimated sample size

#CAREFUL! The r value is not what you want to have in your study, but what you assume the correlation in your study will have.
#also r is the square root of r2, so you can back-calculate what r value you need, because r and r2 do NOT mean the same thing
#Coefficient of correlation is “R” value, r2 = Coefficient of Determination
#r2 shows percentage variation in y which is explained by all the x variables together
#Coefficient of Correlation: is the degree of relationship between two variables say x and y. It can go between -1 and 1.  
#1 indicates that the two variables are moving in unison. -1 means they are opposites.

#Cohen suggests that r values of 0.1, 0.3, and 0.5 represent small, medium, and large effect sizes respectively.
#https://www.statmethods.net/stats/power.html 
