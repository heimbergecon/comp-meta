rm(list = ls()) #clear list

#automatic installation of required packages
packages <- c("xlsx","calibrate","stargazer","sandwich","lmtest","getopt","CausalGAM","ggplot2","reshape2","xts",
              "lattice","gridExtra","gtable","plm","lfe","lmtest","car","tis","foreign","MASS","quantreg","ggrepel",
              "dplyr","stringr","datasets","rio","psych","systemfit","MatchIt","CRTgeeDR","eurostat","plyr","zoo","ggthemes",
              "robumeta","metafor","dplyr","clubSandwich","Hmisc","metafor","pracma","pkgs","broom","sjPlot", "here", "data.table")
ipak(packages)

#load packages
library(xlsx) #Excel-Paket laden
library(calibrate) #Laden des Pakets, das f??r Datenbeschriftung n??tig ist
library (stargazer) #Laden des Pakets, mit dem R-Regressionsoutput in Latex-Tabellen ??bergef??hrt werden kann
library(sandwich)
library(lmtest)
library(getopt)
library(CausalGAM)
library(ggplot2)
library(reshape2)
library(xts)
library(lattice)
library(gridExtra)
library(gtable)
library(plm)
library(lfe)
library(lmtest)
library(car)
library(tis)
library(foreign)
library(MASS)
library(quantreg)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggplot2)
library(datasets)
library(rio)
library(psych)
library(systemfit)
library(foreign)
library(MatchIt)
library(CRTgeeDR)
library(eurostat)
library(plyr)
library(zoo)
library(ggthemes)
library("robumeta")
library("metafor")
library("dplyr")
library(clubSandwich)
library(Hmisc)
library(metafor)
library(pracma)
library(pkgs)
library(broom)
library(sjPlot)
library(here)
library(data.table)

dat <- fread(here("Meta-data-tax-competition-FINAL.csv"))

#calculate the partial correlation coefficient
dat$PartialCorrelationCoefficient <- dat$Tstatistic / (sqrt((dat$Tstatistic^2)+dat$DegreesofFreedom))

dat$PartialCorrelationCoefficient <- dat$PartialCorrelationCoefficient*dat$Transform
dat$Tstatistic <- dat$Tstatistic*dat$Transform

#calculate the standard error of the partial correlation coefficient
dat$StandardErrorPartialCorrelation <- sqrt((1-(dat$PartialCorrelationCoefficient)^2)/dat$DegreesofFreedom)

#Precision
dat$PrecSE <- 1 / dat$StandardErrorPartialCorrelation
#InverseSE
dat$InverseSE <- 1 / dat$StandardError
#PrecSESemiElasticity
dat$PrecSESemiElasticity <- 1 / dat$StandardErrorSemiElasticity
#Variance SemiElasticity
dat$VarianceSemiElasticity <- dat$StandardErrorSemiElasticity^2
#Variance 
dat$Variance <- dat$StandardErrorPartialCorrelation^2
#PrecVariance
dat$PrecVariance <- 1 / dat$Variance
#PrecVariance
dat$PrecVarianceSemiElasticity <- 1 / dat$VarianceSemiElasticity

#Proxy for standard error
dat$proxySE <- 1 / (sqrt(dat$Observations))

dat <- escalc(measure="ZCOR", ri=PartialCorrelationCoefficient, ni=Observations, data=dat)

dat$YearofPublication <- log(dat$YearofPublication)
dat$NumberofCountries <- log(dat$NumberofCountries)

dat$StartYear <- as.numeric(as.integer(dat$StartYear))
dat$EndYear <- as.numeric(as.integer(dat$EndYear))

dat$MeanYearData<- (dat$StartYear+dat$EndYear)/2
dat$MeanYearData <- log(dat$MeanYearData)
dat$InstrumentSE <- 1 / (sqrt(dat$DegreesofFreedom))
dat$InstrumentVariance <- dat$InstrumentSE^2
dat$PrecInstrumentVariance <- 1 / dat$InstrumentVariance

dat_long <- melt(dat, id=1:70)

dat_Statutory <- subset(dat_long, Statutory %in% c('1'))
dat_EATR <- subset(dat_long, EATR %in% c('1'))
dat_EMTR <- subset(dat_long, EMTR %in% c('1'))
dat_AverageTaxRate <- subset(dat_long, AverageTaxRate %in% c('1'))
dat_long_SemiElasticity <- subset(dat_long, SEYN %in% c('1'))
dat_long_SemiElasticity_Statutory <- subset(dat_long_SemiElasticity, Statutory %in% c('1'))
dat_long_SemiElasticity_EATR <- subset(dat_long_SemiElasticity, EATR %in% c('1'))
dat_long_SemiElasticity_EMTR <- subset(dat_long_SemiElasticity, EMTR %in% c('1'))

topbottom <- group_by(dat_long, id) %>%
  mutate(rank = rank(desc(PartialCorrelationCoefficient))) %>%
  filter(rank >=31) %>%
  filter(rank <=573) %>%
  arrange(rank)

#Descriptive statistics (Table 1 for partial correlations; Table A4 for semi-elasticities)

#all estimates
#median
#Partial correlation
median(dat$PartialCorrelationCoefficient, na.rm=TRUE)

#Semi-elasticity
median(dat_long_SemiElasticity$SemiElasticity, na.rm=TRUE)

#Unweighted average
#Partial correlation
reguwa <- lm(PartialCorrelationCoefficient~1, data=dat_long)
summary(reguwa)

#Semi-Elasticity
reguwa_SemiElasticity <- lm(SemiElasticity~1, data=dat_long_SemiElasticity)
summary(reguwa_SemiElasticity)

#WLS
#Partial correlation
regwa <- lm(PartialCorrelationCoefficient~1, data=dat_long, weights=PrecVariance)
summary(regwa)
confint(regwa, level=0.95)

#Semi-elasticity
dat_long_SemiElasticity$PrecisionSESemiElasticity <- 1/dat_long_SemiElasticity$VarianceSemiElasticity
regwa_SemiElasticity <- lm(SemiElasticity~1, data=dat_long_SemiElasticity, weights=PrecisionSESemiElasticity)
summary(regwa_SemiElasticity)
confint(regwa_SemiElasticity, level=0.95)

#top-bottom
#median
#Partial correlation
median(topbottom$PartialCorrelationCoefficient, na.rm=TRUE)

#Unweighted average
#Partial correlation
reguwa <- lm(PartialCorrelationCoefficient~1, data=topbottom)
summary(reguwa)

#WLS
#Partial correlation
regwa <- lm(PartialCorrelationCoefficient~1, data=topbottom, weights=PrecVariance)
summary(regwa)
confint(regwa, level=0.95)

#Statutory
#median
#Partial correlation
median(dat_Statutory$PartialCorrelationCoefficient, na.rm=TRUE)

#Semi-elasticity
median(dat_long_SemiElasticity_Statutory$SemiElasticity, na.rm=TRUE)

#Unweighted average
#Partial correlation
reguwa <- lm(PartialCorrelationCoefficient~1, data=dat_Statutory)
summary(reguwa)

#Semi-elasticity
reguwa_SemiElasticity <- lm(SemiElasticity~1, data=dat_long_SemiElasticity_Statutory)
summary(reguwa_SemiElasticity)

#WLS
#Partial correlation
regwa <- lm(PartialCorrelationCoefficient~1, data=dat_Statutory, weights=PrecVariance)
summary(regwa)
confint(regwa, level=0.95)

#Semi-elasticity
dat_long_SemiElasticity_Statutory$PrecisionSESemiElasticity <- 1/dat_long_SemiElasticity_Statutory$VarianceSemiElasticity
regwa_SemiElasticity <- lm(SemiElasticity~1, data=dat_long_SemiElasticity_Statutory, weights=PrecisionSESemiElasticity)
summary(regwa_SemiElasticity)
confint(regwa_SemiElasticity, level=0.95)

#Semi-elasticity
regwa_re <- rma(SemiElasticity, VarianceSemiElasticity, mods = ~ 1, weights=PrecSE, weighted=FALSE, method="REML", data=dat_long_SemiElasticity_Statutory) 
summary(regwa_re)
confint(regwa_re, level=0.95)

#EATR
#median
#Partial correlation
median(dat_EATR$PartialCorrelationCoefficient, na.rm=TRUE)

#Semi-elasticity
median(dat_long_SemiElasticity_EATR$SemiElasticity, na.rm=TRUE)

#Unweighted average
#Partial correlation
reguwa <- lm(PartialCorrelationCoefficient~1, data=dat_EATR)
summary(reguwa)

#Semi-elasticity
reguwa_SemiElasticity <- lm(SemiElasticity~1, data=dat_long_SemiElasticity_EATR)
summary(reguwa)

#WLS
#Partial correlation
regwa <- lm(PartialCorrelationCoefficient~1, data=dat_EATR, weights=PrecVariance)
summary(regwa)
confint(regwa, level=0.95)

#Semi-elasticity
dat_long_SemiElasticity_EATR$PrecisionSESemiElasticity <- 1/dat_long_SemiElasticity_EATR$VarianceSemiElasticity
regwa_SemiElasticity <- lm(SemiElasticity~1, data=dat_long_SemiElasticity_EATR, weights=PrecisionSESemiElasticity)
summary(regwa_SemiElasticity)
confint(regwa_SemiElasticity, level=0.95)

#EMTR
#median
#Partial correlation
median(dat_EMTR$PartialCorrelationCoefficient, na.rm=TRUE)

#Semi-elasticity
median(dat_long_SemiElasticity_EMTR$SemiElasticity, na.rm=TRUE)

#Unweighted average
#Partial correlation
reguwa <- lm(PartialCorrelationCoefficient~1, data=dat_EMTR)
summary(reguwa)

#Semi-elasticity
reguwa_SemiElasticity <- lm(SemiElasticity~1, data=dat_long_SemiElasticity_EMTR)
summary(reguwa_SemiElasticity)

#WLS
#Partial correlation
regwa <- lm(PartialCorrelationCoefficient~1, data=dat_EMTR, weights=PrecVariance)
summary(regwa)
confint(regwa, level=0.95)

#Semi-elasticity
dat_long_SemiElasticity_EMTR$PrecisionSESemiElasticity <- 1/dat_long_SemiElasticity_EMTR$VarianceSemiElasticity
regwa_SemiElasticity <- lm(SemiElasticity~1, data=dat_long_SemiElasticity_EMTR, weights=PrecisionSESemiElasticity)
summary(regwa_SemiElasticity)
confint(regwa_SemiElasticity, level=0.95)

#ATR
#median
#Partial correlation
median(dat_AverageTaxRate$PartialCorrelationCoefficient, na.rm=TRUE)

#Unweighted average
#Partial correlation
reguwa <- lm(PartialCorrelationCoefficient~1, data=dat_AverageTaxRate)
summary(reguwa)

#WLS
#Partial correlation
regwa <- lm(PartialCorrelationCoefficient~1, data=dat_AverageTaxRate, weights=PrecVariance)
summary(regwa)
confint(regwa, level=0.95)

#minimum, maximum, standard deviation
max(dat_long$PartialCorrelationCoefficient)
min(dat_long$PartialCorrelationCoefficient)
sd(dat_long$PartialCorrelationCoefficient)
mean(dat_long$PartialCorrelationCoefficient)
median(dat_long$PartialCorrelationCoefficient)

#Visualization of publication selection bias
#Figure 2

#Partial correlations
plot_funnel_partialcorrelation <- ggplot(data=dat,
                      aes(x=PartialCorrelationCoefficient, y=PrecSE)) +
  geom_point(aes(colour=factor(DependentVariable)), size=0.5) +
  xlab("Partial correlation coefficient") +
  ylab("Inverse of standard error (precision)") +
  ggtitle("")+
  theme(title=element_text(size=10, face='bold'))+
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank()) +
  geom_vline(xintercept=0, colour="black", linetype=2)+
  geom_vline(xintercept=0.094, colour="black", linetype=1)+
  theme(legend.text = element_text(colour="black", size = 4))+
  theme(axis.text.x=element_text(size=11))+
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.y=element_text(size=11))
plot_funnel_partialcorrelation

#Figure 3
#Semi-elasticities
plot_funnel_semielasticities <- ggplot(data=dat_long_SemiElasticity,
                      aes(x=SemiElasticity, y=PrecSESemiElasticity)) +
  geom_point(aes(colour=factor(DependentVariable)), size=0.5) +
  xlim(-10,20) +
  xlab("Semi-elasticity") +
  ylab("Inverse of standard error (precision)") +
  ggtitle("")+
  theme(title=element_text(size=10, face='bold'))+
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank()) +
  geom_vline(xintercept=0, colour="black", linetype=2)+
  geom_vline(xintercept=0.79, colour="black", linetype=1)+
  theme(legend.text = element_text(colour="black", size = 4))+
  theme(axis.text.x=element_text(size=11))+
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.text.y=element_text(size=11))+
  theme(axis.title.y=element_text(size=11))
plot_funnel_semielasticities

###
#Multivariate meta-regression analysis (all corporate tax variables)

#Full model (all moderator variables)
pubbias_var_est <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + EATR + EMTR + AverageTaxRate + UniformWeight + EconomicTiesWeight + SizeWeight + OpennessWeight + DataNonEuropeancountries + DataMixofCountries + OLS + GMM + OtherEstimator + StackelbergModel + InteractionSpatialLag + CrossSection + MeanYearData + IntraNational + CountryFixedEffects + TimeFixedEffects + TimeTrend + ReviewedJournal + Prior + Crossauthor + OpennessVariableIncluded + OldAgePopulation + Unemployment + PartisanPolitics + Size, weights=PrecVariance, data=dat_long)
summary(pubbias_var_est)

coef_test(pubbias_var_est, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

#general-to-specific
pubbias_var_est_gts <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation + EATR + EMTR + AverageTaxRate + GeographicDistanceWeight + EconomicTiesWeight + SizeWeight + OpennessWeight + DataNonEuropeancountries + DataMixofCountries + OLS + GMM + OtherEstimator + CrossSection + Prior + PartisanPolitics + Size, weights=PrecVariance, data=dat_long)
summary(pubbias_var_est_gts)

coef_test(pubbias_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

#column (3)
pubbias_var_est_gts_IV <- lm(PartialCorrelationCoefficient ~ InstrumentSE + EATR + EMTR + AverageTaxRate + GeographicDistanceWeight + EconomicTiesWeight + SizeWeight + OpennessWeight + DataNonEuropeancountries + DataMixofCountries + OLS + GMM + OtherEstimator + CrossSection + Prior + PartisanPolitics + Size, weights=PrecVariance, data=dat_long)
summary(pubbias_var_est_gts_IV)

coef_test(pubbias_var_est_gts_IV, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")


#column(4)
#Fisher's z
pubbias_var_est_gts_Fishers <- lm(yi ~ StandardErrorPartialCorrelation + EATR + EMTR + AverageTaxRate + GeographicDistanceWeight + EconomicTiesWeight + SizeWeight + OpennessWeight + DataNonEuropeancountries + DataMixofCountries + OLS + GMM + OtherEstimator + CrossSection + Prior + PartisanPolitics + Size, weights=PrecVariance, data=dat_long)
summary(pubbias_var_est_gts_Fishers)

coef_test(pubbias_var_est_gts_Fishers, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

#stargazer table

#WLS without clustered standard errors
ses.WLS.est.nocluster <- list(coeftest(pubbias_var_est_gts)[,2]) #heteroskedasticity-robust standard errors
tvals.WLS.est.nocluster <- list(coeftest(pubbias_var_est_gts)[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.
tvals.WLS.est.nocluster
ses.WLS.est.nocluster

#WLS with standard errors clustered at the study level
#general-to-specific
ses.WLS.est <- list(coef_test(pubbias_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                              cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals.WLS.est <- list(coef_test(pubbias_var_est_gts, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.
tvals.WLS.est
ses.WLS.est

#WLS - full model with standard errors clustered at the study level
#general-to-specific
ses.WLS.est.full <- list(coef_test(pubbias_var_est, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                              cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals.WLS.est.full <- list(coef_test(pubbias_var_est, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.
tvals.WLS.est.full
ses.WLS.est.full

#Fisher's z
ses.WLS_Fishers <- list(coef_test(pubbias_var_est_gts_Fishers, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                              cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals.WLS_Fishers <- list(coef_test(pubbias_var_est_gts_Fishers, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.
tvals.WLS_Fishers
ses.WLS_Fishers

#IV
ses.WLS_IV <- list(coef_test(pubbias_var_est_gts_IV, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                              cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals.WLS_IV <- list(coef_test(pubbias_var_est_gts_IV, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.
tvals.WLS_IV
ses.WLS_IV

#Output for Table 3
stargazer(pubbias_var_est_gts, pubbias_var_est_gts, pubbias_var_est_gts_Fishers, pubbias_var_est_gts_IV, t=list(unlist(tvals.WLS.est.nocluster), unlist(tvals.WLS.est), unlist(tvals.WLS_Fishers), unlist(tvals.WLS_IV)), se=list(unlist(ses.WLS.est.nocluster), unlist(ses.WLS.est), unlist(ses.WLS_Fishers), unlist(ses.WLS_IV)))

#Output for Table A3
stargazer(pubbias_var_est, t=list(unlist(tvals.WLS.est.full)), se=list(unlist(ses.WLS.est.full)))

#descriptive statistics for table 1

mean(dat_long$PartialCorrelationCoefficient)

mean(dat_long$Statutory)
mean(dat_long$EATR)
mean(dat_long$EMTR)
mean(dat_long$AverageTaxRate)

mean(dat_long$UniformWeight)
mean(dat_long$GeographicDistanceWeight)
mean(dat_long$SizeWeight)
mean(dat_long$EconomicTiesWeight)
mean(dat_long$OpennessWeight)

mean(dat_long$CrossSection)
mean(dat_long$IntraNational)
mean(dat_long$MeanYearData)

mean(dat_long$DataEuropeanCountries)
mean(dat_long$DataNonEuropeancountries)
mean(dat_long$DataMixofCountries)

mean(dat_long$StackelbergModel)
mean(dat_long$InteractionSpatialLag)
mean(dat_long$CountryFixedEffects)
mean(dat_long$TimeFixedEffects)
mean(dat_long$TimeTrend)
mean(dat_long$IV2SLS)
mean(dat_long$OLS)
mean(dat_long$GMM)
mean(dat_long$OtherEstimator)

mean(dat_long$StandardErrorPartialCorrelation)
mean(dat_long$ReviewedJournal)
mean(dat_long$Primary)
mean(dat_long$Crossauthor)
mean(dat_long$Prior)

mean(dat_long$OpennessVariableIncluded)
mean(dat_long$OldAgePopulation)
mean(dat_long$Unemployment)
mean(dat_long$PartisanPolitics)
mean(dat_long$GovernmentSpending)
mean(dat_long$Size)

sd(dat_long$PartialCorrelationCoefficient)

sd(dat_long$Statutory)
sd(dat_long$EATR)
sd(dat_long$EMTR)
sd(dat_long$AverageTaxRate)

sd(dat_long$UniformWeight)
sd(dat_long$GeographicDistanceWeight)
sd(dat_long$SizeWeight)
sd(dat_long$EconomicTiesWeight)
sd(dat_long$OpennessWeight)

sd(dat_long$CrossSection)
sd(dat_long$IntraNational)
sd(dat_long$MeanYearData)

sd(dat_long$DataEuropeanCountries)
sd(dat_long$DataNonEuropeancountries)
sd(dat_long$DataMixofCountries)

sd(dat_long$StackelbergModel)
sd(dat_long$InteractionSpatialLag)
sd(dat_long$CountryFixedEffects)
sd(dat_long$TimeFixedEffects)
sd(dat_long$TimeTrend)
sd(dat_long$IV2SLS)
sd(dat_long$OLS)
sd(dat_long$GMM)
sd(dat_long$OtherEstimator)

sd(dat_long$StandardErrorPartialCorrelation)
sd(dat_long$ReviewedJournal)
sd(dat_long$Primary)
sd(dat_long$Crossauthor)
sd(dat_long$Prior)

sd(dat_long$OpennessVariableIncluded)
sd(dat_long$OldAgePopulation)
sd(dat_long$Unemployment)
sd(dat_long$PartisanPolitics)
sd(dat_long$GovernmentSpending)
sd(dat_long$Size)

#
#Results on publication selection bias (supplementary appendix, table A2)

#all corporate tax variables
#average precision-weighted partial correlation
average_effect <- lm(PartialCorrelationCoefficient ~1, weights=PrecVariance, data=dat_long)
summary(average_effect)
coef_test(average_effect, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

#column (1)
#FAT-PET test
pubbias_1 <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_long)
summary(pubbias_1)

coef_test(pubbias_1, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_long$id, test = "naive-t")

dat_long$PESEE <- dat$StandardErrorPartialCorrelation^2

#column (2): Statutory tax rates
pubbias_3_Statutory <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_Statutory)
summary(pubbias_3_Statutory)

coef_test(pubbias_3_Statutory, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_Statutory$id, test = "naive-t")

#column (3): EATR
pubbias_3_EATR <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_EATR)
summary(pubbias_3_EATR)

coef_test(pubbias_3_EATR, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_EATR$id, test = "naive-t")

#column (4): EMTR
pubbias_3_EMTR <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_EMTR)
summary(pubbias_3_EMTR)

coef_test(pubbias_3_EMTR, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_EMTR$id, test = "naive-t")

#column (5): AverageTaxRate
pubbias_3_AverageTaxRate <- lm(PartialCorrelationCoefficient ~ StandardErrorPartialCorrelation, weights=PrecVariance, data=dat_AverageTaxRate)
summary(pubbias_3_AverageTaxRate)

coef_test(pubbias_3_AverageTaxRate, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
          cluster = dat_AverageTaxRate$id, test = "naive-t")

#correct t-values and standard errors for stargazer table
#all corporate tax variable
ses_pubbias_1 <- list(coef_test(pubbias_1, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                cluster = dat_long$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_1 <- list(coef_test(pubbias_1, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                  cluster = dat_long$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#Statutory
ses_pubbias_3_Statutory <- list(coef_test(pubbias_3_Statutory, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                          cluster = dat_Statutory$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_3_Statutory <- list(coef_test(pubbias_3_Statutory, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                            cluster = dat_Statutory$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


#EATR
ses_pubbias_3_EATR <- list(coef_test(pubbias_3_EATR, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                     cluster = dat_EATR$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_3_EATR <- list(coef_test(pubbias_3_EATR, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                       cluster = dat_EATR$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#EMTR
ses_pubbias_3_EMTR <- list(coef_test(pubbias_3_EMTR, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                     cluster = dat_EMTR$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_3_EMTR <- list(coef_test(pubbias_3_EMTR, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                       cluster = dat_EMTR$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.

#AverageTaxRate
ses_pubbias_3_AverageTaxRate <- list(coef_test(pubbias_3_AverageTaxRate, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                               cluster = dat_AverageTaxRate$id, test = "naive-t")[,2]) #heteroskedasticity-robust standard errors
tvals_pubbias_3_AverageTaxRate <- list(coef_test(pubbias_3_AverageTaxRate, vcov = "CR1", #CR1 refers to small sample correction; see package documentation clubSandwich
                                                 cluster = dat_AverageTaxRate$id, test = "naive-t")[,3]) #heteroskedasticity-robust t-values "group" ("time") accounts for serial (cross-sectional) correlation.


stargazer(pubbias_1, pubbias_3_Statutory, pubbias_3_EATR, pubbias_3_EMTR, pubbias_3_AverageTaxRate, t=list(unlist(tvals_pubbias_1), unlist(tvals_pubbias_3_Statutory), unlist(tvals_pubbias_3_EATR), unlist(tvals_pubbias_3_EMTR), unlist(tvals_pubbias_3_AverageTaxRate)), se=list(unlist(ses_pubbias_1), unlist(ses_pubbias_3_Statutory), unlist(ses_pubbias_3_EATR), unlist(ses_pubbias_3_EMTR), unlist(ses_pubbias_3_AverageTaxRate)))

