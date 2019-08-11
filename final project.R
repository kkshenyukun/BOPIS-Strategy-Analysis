---
title: "Final Project"
author: "Yukun Shen"
output: word_document

rm(list = ls())
setwd("/Users/yukunshen/Desktop/BOPS")
library(stargazer)
library(gdata)
library(ggplot2)
library(psych) 
library(ggeffects)
library(QuantPsyc)
library(usdm)
library(lmtest)
library(multiwayvcov)
library(sandwich)
library(foreign)
library(AER)
library(aod)
library(Rcpp)
library(mfx)
library(nnet)
library(reshape2)
library(VIF)
library(MASS)
library(readstata13)
library(ivprobit)
library(Rmisc)
#==========================================================
## Q1: (OLS REGRESSION )
#==========================================================

data1 = read.dta13("/Users/yukunshen/Desktop/BOPS/consumer level data.dta")
data2 = read.dta13("/Users/yukunshen/Desktop/BOPS/online daily prod_cat sales-returns data.dta")
data3 = read.dta13("/Users/yukunshen/Desktop/BOPS/online daily sales-returns data.dta")
data4 = read.dta13("/Users/yukunshen/Desktop/BOPS/transaction level data.dta")
#choose data3 as our dataset
stargazer(data3, type="text", median=TRUE, iqr=TRUE,digits=1, title="Descriptive Statistics")  

ggplot(data3, aes(x=salesvalue)) + geom_histogram(colour="green")
ggplot(data3, aes(x=log(salesvalue))) + geom_histogram(colour="green")#use log transformed salesvalue since the distribution looks more normal
summary(data3)
mydata1 <- subset(data3,day<786)#generate a new data frame that exclude days after 785(when store 5998 start to use bops)
mydata1$group<-ifelse(mydata1$store_number==5998,0,1)#seperate stores by creating a dummy variable
mydata1$time<-ifelse(mydata1$day<366,0,1)#creating a time dummy variable
#make all NA values in our data to mean
mydata1$avg_femalemean<-ifelse(is.na(mydata1$avg_female),mean(mydata1$avg_female,na.rm=TRUE),mydata1$avg_female)
mydata1$avg_incomemean<-ifelse(is.na(mydata1$avg_income),mean(mydata1$avg_income,na.rm=TRUE),mydata1$avg_income)
mydata1$avg_homeownermean<-ifelse(is.na(mydata1$avg_homeowner),mean(mydata1$avg_homeowner,na.rm=TRUE),mydata1$avg_homeowner)
mydata1$avg_residencymean<-ifelse(is.na(mydata1$avg_residency),mean(mydata1$avg_residency,na.rm=TRUE),mydata1$avg_residency)
mydata1$avg_childownermean<-ifelse(is.na(mydata1$avg_childowner),mean(mydata1$avg_childowner,na.rm=TRUE),mydata1$avg_childowner)
mydata1$avg_agemean<-ifelse(is.na(mydata1$avg_age),mean(mydata1$avg_age,na.rm=TRUE),mydata1$avg_age)
#check multicollinearity
df <- mydata1[c("time","group","avg_femalemean","avg_incomemean","avg_homeownermean","avg_residencymean","avg_agemean")]
cor(df) 
vifcor(df)
#build our initial model 
m1 <- lm(log(salesvalue+1)~time*group,data=mydata1)
stargazer(m1, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
is.factor(mydata1$month_dummy)
mydata1$month_dummy <- as.factor(mydata1$month_dummy)
m2 <- lm(log(salesvalue+1)~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=mydata1)
stargazer(m2, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 

#check which model is better
anova(m1, m2, test="Chisq")   #significant P-value indicates our second model is better

# Check for heteroscedasticity
gqtest(m2)
bptest(m2)
HWrobstder <- sqrt(diag(vcovHC(m2, type="HC1")))
stargazer(m2,m2,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))

#==========================================================
## Q1:POISSON REGRESSION & Negative Binomial
#==========================================================
#we can't use linear probability model because our dependent variable is a count variable.By running an linear probability model we would have negative sales quantity.

#build our initial poisson model
poisson1 <- glm(salesquantity~time*group,family="poisson",data=mydata1)
stargazer(poisson1,  
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))
poisson1a <- glm(salesquantity~1, data=mydata1, family="poisson") #run a logit on null model 
lrtest(poisson1, poisson1a)#after running a likelihood test, we found a significant P-Value. Thus, we conclude our data does not fit the model
#add control variables
poisson2 <- glm(salesquantity~time*group+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean,family="poisson",data=mydata1)
stargazer(poisson1,  poisson2,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson2, poisson1a)#after running a likelihood test, we found a significant P-Value. Thus, we conclude our data does not fit the model

#NEGATIVE BINOMIAL MODEL
negbin1 <- glm.nb(salesquantity~time*group,data=mydata1)

stargazer(negbin1,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin2 <- glm.nb(salesquantity~time*group+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean,data=mydata1)

stargazer(negbin1,  negbin2,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))
 
negbin1a <- glm.nb(log(salesquantity) ~ 1, data = mydata1) #build a null model

lrtest(negbin2, negbin1a) #significant P-value means model fit
#choose model
lrtest(poisson2,negbin2)#siginificant p-value means negative binomial model is more appropriate

#check for heterokadasticity
gqtest(negbin2) 
bptest(negbin2)
HWrobstder <- sqrt(diag(vcovHC(negbin2, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin2,negbin2,
          se=list(NULL, HWrobstder),
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  
#==========================================================
##Q2 (OLS)
#==========================================================
#use the same dataset as Q1
ggplot(data3, aes(x=returnvalue)) + geom_histogram(colour="green")
ggplot(data3, aes(x=log(returnvalue))) + geom_histogram(colour="green") #use log transformed return value since the distribution looks more normal

#build initial model
m21 <- lm(log(returnvalue+1)~time*group,data=mydata1)
stargazer(m21, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
m22 <- lm(log(returnvalue+1)~time*group+avg_agemean+avg_femalemean+avg_incomemean+log(salesvalue+1)+month_dummy+avg_childownermean+avg_homeownermean,data=mydata1)
stargazer(m21, m22,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))     
#check which model is better
anova(m21,m22)
# Check for heteroscedasticity
gqtest(m22)
bptest(m22)
HWrobstde <- sqrt(diag(vcovHC(m22, type="HC1")))
stargazer(m22,m22,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))
#==========================================================
##Q2 (POISSON REGRESSION & Negative Binomial)
#==========================================================
#we can't use linear probability model because our dependent variable is a count variable.By running an linear probability model we would have negative return quantity.
#build initial poisson model
poisson21 <- glm(returnquantity~time*group,family="poisson",data=mydata1)
stargazer(poisson21,  
          apply.coef = exp, t.auto=F, p.auto = F,
           title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson21a <- glm(returnquantity~1, data=mydata1, family="poisson") #run a logit on null model 
lrtest(poisson1, poisson1a)#siginificant p-value indicates no fit
#add control variables
poisson22 <- glm(returnquantity~time*group+avg_agemean+avg_incomemean+avg_femalemean+salesquantity+month_dummy+avg_childownermean,family="poisson",data=mydata1)
stargazer(poisson21, poisson22,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
anova(poisson21,poisson22,test="Chisq")#indicates poisson 22 is better
lrtest(poisson22, poisson21a)#significant p-value indicates no fit

#NEGATIVE BINOMIAL MODEL
negbin21 <- glm.nb(returnquantity~time*group,data=mydata1)

stargazer(negbin21,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))

negbin22 <- glm.nb(returnquantity~time*group+avg_femalemean+avg_incomemean+salesquantity+month_dummy+avg_agemean+avg_childownermean,data=mydata1)

stargazer(negbin21,  negbin22,
          apply.coef = exp, t.auto=F, p.auto = F,
           title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#generate an empty model
negbin21a <- glm.nb(returnquantity ~ 1, data = mydata1) 
lrtest(negbin22, negbin21a)  #P-value significant means model fit
#choose which model is better
lrtest(poisson22,negbin22)  #significant P-value indicates negtive binomial model is better

#check for heterokadasticity
gqtest(negbin22) 
bptest(negbin22)
HWrobstder <- sqrt(diag(vcovHC(negbin22, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin22,negbin22,
          se=list(NULL, HWrobstder),
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  


#==========================================================
##Q3 (OLS)
#==========================================================
#use data1, which represents consumer level data
summary(data1)
ggplot(data1, aes(x=salesvalue)) + geom_histogram(colour="black")
ggplot(data1, aes(x=log(salesvalue))) + geom_histogram(colour="black")
#decide to use log transformed salesvalue because the distribution looks more normal

#change all NA values into median 
data1$est_income_codemedian<-ifelse(is.na(data1$est_income_code),median(data1$est_income_code,na.rm=TRUE),data1$est_income_code)
data1$age_bandmedian<-ifelse(is.na(data1$age_band),median(data1$age_band,na.rm=TRUE),data1$age_band)

#change child & homeowner_code into dummy variable
data1$childdummy <- ifelse(data1$child=="Y",1,0)
data1$homeowner_codedummy <- ifelse(data1$homeowner_code=="O",1,0)

#generate initial model
m31 <- lm(log(salesvalue+1)~bops_user*bops_in_effect,data=data1)
stargazer(m31, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
m32 <- lm(log(salesvalue+1)~bops_user*bops_in_effect+est_income_codemedian+age_bandmedian+childdummy+homeowner_codedummy+purchase_time_period,data=data1)
stargazer(m31, m32,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 
# Check for heteroscedasticity
gqtest(m32)
bptest(m32)
HWrobstder <- sqrt(diag(vcovHC(m32, type="HC1")))
stargazer(m32,m32,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))

#==========================================================
##Q3 Poisson & Negative Binomial
#==========================================================
#build initial poisson model
poisson31 <- glm(salesquantity~bops_user*bops_in_effect,family="poisson",data=data1)
stargazer(poisson31,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson31a <- glm(salesquantity~1, data=data1, family="poisson") #run a logit on null model 
lrtest(poisson31, poisson31a)#signigicant p-value indicates no fit
#add control variables
poisson32 <- glm(salesquantity~bops_user*bops_in_effect+est_income_codemedian+age_bandmedian+childdummy+purchase_time_period+homeowner_codedummy,family="poisson",data=data1)
stargazer(poisson31,  poisson32,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson32, poisson31a)#significant p-value indicates no fit

#NEGATIVE BINOMIAL MODEL
negbin31 <- glm.nb(salesquantity~bops_user*bops_in_effect,data=data1)

stargazer(negbin31,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin32 <- glm.nb(salesquantity~bops_user*bops_in_effect+est_income_codemedian+age_bandmedian+childdummy+purchase_time_period+homeowner_codedummy,data=data1)

stargazer(negbin31,  negbin32,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))

negbin31a <- glm.nb(log(salesquantity) ~ 1, data = data1) #generate a null model for comparison

lrtest(negbin32, negbin31a) #significant P-value means model fit
#choose model
lrtest(poisson32,negbin32) #significant P-Value indicates negative binomial is better

# Check for heteroscedasticity
gqtest(negbin32)
bptest(negbin32)
HWrobstder <- sqrt(diag(vcovHC(negbin32, type="HC1")))
stargazer(negbin32,negbin32,
          apply.coef = exp, t.auto=F, p.auto = F,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#==========================================================
##Q4
#==========================================================
#use data 4 which represents transaction level data
mydata4 <- data4[!is.na(data4$bops), ] #exclude all NA variables in bops
#transform specific variables into dummy variable and factory variable
mydata4$childdummy <- ifelse(mydata4$child=="Y",1,0)
mydata4$homeowner_codedummy <- ifelse(mydata4$homeowner_code=="O",1,0)
mydata4$est_income_codemedian<-ifelse(is.na(mydata4$est_income_code),median(mydata4$est_income_code,na.rm=TRUE),mydata4$est_income_code)
mydata4$age_bandmedian<-ifelse(is.na(mydata4$age_band),median(mydata4$age_band,na.rm=TRUE),mydata4$age_band)
mydata4$length_of_residencemedian<-ifelse(is.na(mydata4$length_of_residence),median(mydata4$length_of_residence,na.rm=TRUE),mydata4$length_of_residence)
#excluding NA/NULL values
mydata41 <- mydata4[!is.na(mydata4$female),]
mydata42<-mydata41[!mydata41$child=="", ]
mydata42$childdummy <- ifelse(mydata42$child=="Y",1,0)
mydata42$homeowner_codedummy <- ifelse(mydata42$homeowner_code=="O",1,0)

ggplot(mydata42, aes(x=price)) + geom_histogram(colour="black")
ggplot(mydata42, aes(x=log(price))) + geom_histogram(colour="black")#log transformed price looks more normal
#log transformed price and change several variables into factor variable
mydata42$logprice <- log(mydata42$price+1)
mydata42$storefactor <- as.factor(mydata42$store_number)
mydata42$monthfactor <- as.factor(mydata42$month_dummy)
mydata42$yearfactor <- as.factor(mydata42$year)
stargazer(mydata42,type="text",median=TRUE,iqr=TRUE,digits=1,title = "Descriptive Statistics")
#build initial OLS model
m41 <- lm(return~bops+price+age_band+est_income_code+female,data=mydata42)
stargazer(m41, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 
#using new dataset with replaced median variable
m41a <- lm(return~bops+logprice+age_bandmedian+est_income_codemedian+female+storefactor+monthfactor+yearfactor+childdummy,data=mydata42)
stargazer(m41a, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))

#LOGIT MODEL
sum(mydata42$return==0)
sum(mydata42$return==1) # We have 99452 observations with Return=1 and 887466 observations with Return=0. Considering that we will estimate 29 parameters, we satisfy the minimum 10:1 ratio requirement
#build initial logit model
logit41 <- glm(return~bops+logprice+age_bandmedian+est_income_codemedian+female+storefactor+monthfactor+yearfactor+childdummy,data=mydata42,family="binomial")
stargazer(logit41, m41a,
          title="Regression Results", type="text", 
          column.labels=c("Logit-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))
stargazer(logit41, m41a,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("OddsRatios"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001))
logit1a <- glm(return~1, data=mydata42, family="binomial") # run a logit on null model 
lrtest(logit41, logit1a)#comparing our model with null model

#generate marginal effects
a <- logitmfx(formula=return~bops+logprice+age_bandmedian+est_income_codemedian+female+storefactor+monthfactor+yearfactor+childdummy,data=mydata42) 
marginaleffects <- a$mfxest[,1]
marg.std.err <- a$mfxest[,2]

stargazer(logit41,
          omit=c("Constant"),
          coef = list(marginaleffects), se = list(marg.std.err),
          title="Regression Results", type="text", 
          column.labels=c("Marginal Effects"),
          df=FALSE, digits=5, star.cutoffs = c(0.05,0.01,0.001))
#after generate marginal effects, our beta for bops turned out to be 0.016,which is very close to our OLS beta(0.015). Thus we can use IV/2SLS model

#build IV/2SLS model
df4 <- mydata42[c("bops","logprice","age_bandmedian","est_income_codemedian","female","childdummy","length_of_residencemedian","homeowner_codedummy")]
cor(df4)
vifcor(df4)
iv41 <- ivreg(return~bops+logprice+age_bandmedian+est_income_codemedian+female+storefactor+monthfactor+yearfactor|length_of_residencemedian+homeowner_codedummy+logprice+age_bandmedian+est_income_codemedian+female+storefactor+monthfactor+yearfactor,data=mydata42)
summary(iv41)
stargazer(iv41,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
summary(iv41,diagnostics = TRUE)#p-value for Sargan statistic is significant means we can not use our instrumental variable

#generate a dataframe contain only endogenous variable and instrument variable
df41 <- mydata42[c("return","length_of_residencemedian","homeowner_codedummy")]
cor(df41)#drop homeowner_codedummy because of the higher correlation
iv42 <- ivreg(return~bops+logprice+age_bandmedian+est_income_codemedian+female+storefactor+monthfactor+yearfactor|length_of_residencemedian+logprice+age_bandmedian+est_income_codemedian+female+storefactor+monthfactor+yearfactor,data=mydata42)
summary(iv42)
stargazer(iv42,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 
summary(iv42,diagnostics = TRUE)#because this is a just identified model we do not have sargan values. Thus, we have to use expertise knowldge to pass this test. Then significant p-value of Hausman test indicates model fit

# Check for heteroscedasticity
gqtest(iv42)
bptest(iv42)
HWrobstder <- sqrt(diag(vcovHC(iv42, type="HC1")))
stargazer(iv42,iv42,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#==========================================================
##Q5 SALESVALUE
#==========================================================
#use data3 which represents online daily prod_cat sales-returns data
#SALESVALUE
stargazer(data2, type="text", median=TRUE, iqr=TRUE,digits=1, title="Descriptive Statistics")  
ggplot(data2, aes(x=salesvalue)) + geom_histogram(colour="green")
ggplot(data2, aes(x=log(salesvalue))) + geom_histogram(colour="green")#use log transformed salesvalue since the distribution looks more normal
summary(data2)
mydata5 <- subset(data2,day<786)#generate a new data frame that exclude days after 785
mydata5$group<-ifelse(mydata5$store_number==5998,0,1)#seperate stores by creating a dummy variable
mydata5$time<-ifelse(mydata5$day<366,0,1)#creating a time dummy variable
mydata5$avg_femalemean<-ifelse(is.na(mydata5$avg_female),mean(mydata5$avg_female,na.rm=TRUE),mydata5$avg_female)
mydata5$avg_incomemean<-ifelse(is.na(mydata5$avg_income),mean(mydata5$avg_income,na.rm=TRUE),mydata5$avg_income)
mydata5$avg_homeownermean<-ifelse(is.na(mydata5$avg_homeowner),mean(mydata5$avg_homeowner,na.rm=TRUE),mydata5$avg_homeowner)
mydata5$avg_residencymean<-ifelse(is.na(mydata5$avg_residency),mean(mydata5$avg_residency,na.rm=TRUE),mydata5$avg_residency)
mydata5$avg_childownermean<-ifelse(is.na(mydata5$avg_childowner),mean(mydata5$avg_childowner,na.rm=TRUE),mydata5$avg_childowner)
mydata5$avg_agemean<-ifelse(is.na(mydata5$avg_age),mean(mydata5$avg_age,na.rm=TRUE),mydata5$avg_age)
df5<-mydata5[c("time","avg_femalemean","avg_incomemean","avg_homeownermean","avg_residencymean","avg_childownermean","avg_agemean","group","returnvalue","returnquantity","product_category")]
cor(df5) 
vifcor(df5)
m51 <- lm(log(salesvalue+1)~time*group,data=mydata5)
stargazer(m51, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
is.factor(mydata5$month_dummy)
mydata5$month_dummy <- as.factor(mydata5$month_dummy)#change month into factor variable
is.factor(mydata5$product_category)
mydata5$product_category<-as.factor(mydata5$product_category)
m52 <- lm(log(salesvalue+1)~time*group+product_category+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=mydata5)
stargazer(m51,m52, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 

anova(m51,m52)#significant p-value idicates our second model is better
# Check for heteroscedasticity
gqtest(m52)
bptest(m52)
HWrobstder <- sqrt(diag(vcovHC(m52, type="HC1")))
stargazer(m52,m52,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#=========================================================
##Q5 RETURNVALUE
#==========================================================
ggplot(mydata5, aes(x=returnvalue)) + geom_histogram(colour="green")
ggplot(mydata5, aes(x=log(returnvalue))) + geom_histogram(colour="green")#use log transformed returnvalue since the distribution looks more normal

#Build initial model
m51r <- lm(log(returnvalue+1)~time*group,data=mydata5)
stargazer(m51r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 

#add control variables
m52r <- lm(log(returnvalue+1)~time*group+log(salesvalue+1)+product_category+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=mydata5)
stargazer(m51r,m52r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
anova(m51r,m52r)#significant p-value idicates our second model is better
#Check for heterokedasticity
gqtest(m52r)
bptest(m52r)
HWrobstder <- sqrt(diag(vcovHC(m52r, type="HC1")))
stargazer(m52r,m52r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#==========================================================
##Q5 SALESQUANTITY
#==========================================================
poisson51 <- glm(salesquantity~time*group,family="poisson",data=mydata5)
stargazer(poisson51,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson51a <- glm(salesquantity~1, data=mydata5, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson51, poisson51a)#siginificant p-value indicates no fit
#add control variables
poisson52 <- glm(salesquantity~time*group+product_category+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean,family="poisson",data=mydata5)
stargazer(poisson51,  poisson52,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson52, poisson51a)

#NEGATIVE BINOMIAL MODEL
negbin51 <- glm.nb(salesquantity~time*group,data=mydata5)

stargazer(negbin51,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin52 <- glm.nb(salesquantity~time*group+product_category+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean,data=mydata5)

stargazer(negbin51,  negbin52,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))

negbin51a <- glm.nb(salesquantity ~ 1, data = mydata5) #genarate an null model

lrtest(negbin52, negbin51a) #significant P-value means model fit
#choose model
lrtest(poisson52,negbin52)

#check for heterokadasticity
gqtest(negbin52) 
bptest(negbin52)
HWrobstder <- sqrt(diag(vcovHC(negbin52, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin52,negbin52,
          se=list(NULL, HWrobstder),
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  # 

#==========================================================
##Q5 RETURNQUANTITY
#==========================================================

poisson51r <- glm(returnquantity~time*group,family="poisson",data=mydata5)
stargazer(poisson51r,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson51ar <- glm(returnquantity~1, data=mydata5, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson51r, poisson51ar)#significant p-value indicates no fit
#add control variables
poisson52r <- glm(returnquantity~time*group+product_category+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean+salesquantity,family="poisson",data=mydata5)
stargazer(poisson51r,  poisson52r,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson52r, poisson51ar)

#NEGATIVE BINOMIAL MODEL
negbin51r <- glm.nb(returnquantity~time*group,data=mydata5)

stargazer(negbin51r,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin52r <- glm.nb(returnquantity~time*group+product_category+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean+salesquantity,data=mydata5)

stargazer(negbin51r,  negbin52r,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))

negbin51ar <- glm.nb(returnquantity ~ 1, data = mydata5) #generate a null model

lrtest(negbin52r, negbin51ar) #significant P-value means model fit
#choose model
lrtest(poisson52r,negbin52r) #significant p-value means we should use negtiva binomial model

#check for heterokadasticity
gqtest(negbin52r) 
bptest(negbin52r)
HWrobstder <- sqrt(diag(vcovHC(negbin52r, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin52r,negbin52r,
          se=list(NULL, HWrobstder),
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  # 


#==========================================================
##Q6 SALEVALUE
#==========================================================
#use same dataset as last question
mydata5$logsale <- log(mydata5$salesvalue+1)#log transformed salesvalue
mydata5$logreturnvalue <-log(mydata5$returnvalue+1)#use log transformed return value
#build initial model
m61 <- lm(logsale~time*group*product_category,data=mydata5)
stargazer(m61, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
is.factor(mydata5$product_category)
mydata5$product_category <- as.factor(mydata5$product_category)
m62 <- lm(logsale~time*group*product_category+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=mydata5)
stargazer(m61,m62, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 

meffectssalesv <- ggpredict(m62, terms=c("time", "group","product_category")) # generates a tidy data frame  

ggplot(meffectssalesv,aes(x, predicted, colour=group)) + geom_line(size=1.3) + facet_wrap(~facet)
  xlab("bops_in_effect") + ylab("Salesvalue") +
  labs(colour="Payment") + 
  scale_colour_discrete(labels=c("Product Category1", "Product Category")) +
  scale_x_continuous(breaks=c(0,1), labels=c("no bops", "bops")) +
  theme(axis.title.x=element_blank())


# Check for heteroscedasticity
gqtest(m62)
bptest(m62)
HWrobstder <- sqrt(diag(vcovHC(m62, type="HC1")))
stargazer(m62,m62,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#build model for each product category who was affected by bops

#subset product category data
meffects1 <- subset(mydata5,product_category=="1")
meffects2 <- subset(mydata5,product_category=="2")
meffects3 <- subset(mydata5,product_category=="3")
meffects4 <- subset(mydata5,product_category=="4")
meffects5 <- subset(mydata5,product_category=="5")
meffects6 <- subset(mydata5,product_category=="6")
meffects7 <- subset(mydata5,product_category=="7")
meffects8 <- subset(mydata5,product_category=="8")
meffects9 <- subset(mydata5,product_category=="9")
meffects10 <- subset(mydata5,product_category=="10")
meffects11 <- subset(mydata5,product_category=="11")
meffects12 <- subset(mydata5,product_category=="12")
meffects13 <- subset(mydata5,product_category=="13")
meffects14 <- subset(mydata5,product_category=="14")
meffects15 <- subset(mydata5,product_category=="15")
meffects17 <- subset(mydata5,product_category=="17")
meffects20 <- subset(mydata5,product_category=="20")
meffects21 <- subset(mydata5,product_category=="21")

#bops effect on product category 1
me61 <-lm(logsale~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=meffects1)
stargazer(me61, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 

# Check for heteroscedasticity
gqtest(me61)
bptest(me61)
HWrobstder <- sqrt(diag(vcovHC(me61, type="HC1")))
stargazer(me61,me61,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#bops effect on product category 4
me64 <-lm(logsale~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=meffects4)
stargazer(me64, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 

# Check for heteroscedasticity
gqtest(me64)
bptest(me64)
HWrobstder <- sqrt(diag(vcovHC(me64, type="HC1")))
stargazer(me64,me64,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#Bops effect on product category 2
meffects2 <- subset(mydata5,product_category=="2")
me62 <-lm(logsale~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=meffects2)
stargazer(me62, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 

# Check for heteroscedasticity
gqtest(me62)
bptest(me62)
HWrobstder <- sqrt(diag(vcovHC(me64, type="HC1")))
stargazer(me62,me62,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#bops effect on product category 21
meffects21 <- subset(mydata5,product_category=="21")
me621 <-lm(logsale~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=meffects21)
stargazer(me621, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 

# Check for heteroscedasticity
gqtest(me621)
bptest(me621)
HWrobstder <- sqrt(diag(vcovHC(me621, type="HC1")))
stargazer(me621,me621,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#bops effect on product category 14
me614 <-lm(logsale~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean,data=meffects14)
stargazer(me614, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 

# Check for heteroscedasticity
gqtest(me614)
bptest(me614)
stargazer(me614,me614,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#==========================================================
##Q6 RETURNVALUE
#==========================================================
mydata5$logreturnvalue <-log(mydata5$returnvalue+1)#use log transformed return value
m61r <- lm(logreturnvalue~time*group*product_category,data=mydata5)
stargazer(m61r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 
#add control variables
m62r <- lm(logreturnvalue~time*group*product_category+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean+logsale,data=mydata5)
stargazer(m61r,m62r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001)) 

meffectsreturnv <- ggpredict(m62r, terms=c("time", "group","product_category")) # generates a tidy data frame  

ggplot(meffectsreturnv,aes(x, predicted, colour=group)) + geom_line(size=1.3) + facet_wrap(~facet)
xlab("bops_in_effect") + ylab("returnvalue") +
  labs(colour="Payment") + 
  scale_colour_discrete(labels=c("Product Category1", "Product Category")) +
  scale_x_continuous(breaks=c(0,1), labels=c("no bops", "bops")) +
  theme(axis.title.x=element_blank())

# Check for heteroscedasticity
gqtest(m62r)
bptest(m62r)
HWrobstder <- sqrt(diag(vcovHC(m62r, type="HC1")))
stargazer(m62r,m62r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#bops effect on returnvalue product category 4
me64r <-lm(logreturnvalue~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean+logsale,data=meffects4)
stargazer(me64r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 

# Check for heteroscedasticity
gqtest(me64r)
bptest(me64r)
stargazer(me64r,me64r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#bops effect on returnvalue product category 13
me613r <-lm(logreturnvalue~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean+logsale,data=meffects13)
stargazer(me613r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 
# Check for heteroscedasticity
gqtest(me613r)
bptest(me613r)
stargazer(me613r,me613r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#bops effect on returnvalue product category 20
me620r <-lm(logreturnvalue~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean+logsale,data=meffects20)
stargazer(me620r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 

# Check for heteroscedasticity
gqtest(me620r)
bptest(me620r)
stargazer(me620r,me620r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))

#bops effect on returnvalue product category 21
me621r <-lm(logreturnvalue~time*group+avg_agemean+avg_femalemean+avg_incomemean+avg_homeownermean+month_dummy+avg_childownermean+logsale,data=meffects21)
stargazer(me621r, 
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=4, star.cutoffs = c(0.05,0.01,0.001)) 

# Check for heteroscedasticity
gqtest(me621r)
bptest(me621r)
stargazer(me621r,me621r,
          se=list(NULL, HWrobstder),
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))


#==========================================================
##Q6 SALESQUANTITY
#==========================================================
poisson61 <- glm(salesquantity~time*group*product_category,family="poisson",data=mydata5)
stargazer(poisson61,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson61a <- glm(salesquantity~1, data=mydata5, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson61, poisson61a)#significant p-value mean no fit
#add control variables
poisson62 <- glm(salesquantity~time*group*product_category+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean,family="poisson",data=mydata5)
stargazer(poisson61,  poisson62,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson62, poisson61a)#significant p-value mean no fit

#negative binomial model
negbin61 <- glm.nb(salesquantity~time*group*product_category,data=mydata5)

stargazer(negbin61,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin62 <- glm.nb(salesquantity~time*group*product_category+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean,data=mydata5)

stargazer(negbin61,  negbin62,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))

meffectssalesq <- ggpredict(negbin62, terms=c("time", "group","product_category")) # generates a tidy data frame  

ggplot(meffectssalesq,aes(x, predicted, colour=group)) + geom_line(size=1.3) + facet_wrap(~facet)
xlab("bops_in_effect") + ylab("salesquantity") +
  labs(colour="Payment") + 
  scale_colour_discrete(labels=c("Product Category1", "Product Category")) +
  scale_x_continuous(breaks=c(0,1), labels=c("no bops", "bops")) +
  theme(axis.title.x=element_blank())


#check for heterokadasticity
gqtest(negbin62) 
bptest(negbin62)
HWrobstder <- sqrt(diag(vcovHC(negbin62, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin62,negbin62,
          se=list(NULL, HWrobstder),
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  # 
negbin61a <- glm.nb(salesquantity ~ 1, data = mydata5) #Genarate an empty model

lrtest(negbin62, negbin61a) #significant P-value means model fit
#choose model
lrtest(poisson62,negbin62)


#==========================================================
##Q6 RETURNQUANTITY
#==========================================================

poisson61r <- glm(returnquantity~time*group*product_category,family="poisson",data=mydata5)
stargazer(poisson61r,  
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
poisson61ar <- glm(returnquantity~1, data=mydata5, family="poisson") # This is the command to run a logit on null model 
lrtest(poisson61r, poisson61ar)
#add control variables
poisson62r <- glm(returnquantity~time*group*product_category+avg_agemean+avg_femalemean+avg_incomemean+month_dummy+avg_childownermean+avg_homeownermean+salesquantity,family="poisson",data=mydata5)
stargazer(poisson61r,  poisson62r,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(poisson62r, poisson61ar)#significant p-value means no fit

#negative binomial model
negbin61r <- glm.nb(returnquantity~time*group*product_category,data=mydata5)

stargazer(negbin61r,  
          apply.coef = exp, t.auto=F, p.auto = F,
           title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
#add controp varaiables
negbin62r <- glm.nb(returnquantity~time*group*product_category+avg_femalemean+avg_incomemean+avg_agemean+month_dummy+avg_childownermean+avg_homeownermean+salesquantity,data=mydata5)

stargazer(negbin61r,  negbin62r,
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Model-1"),
          df=FALSE, digits=2, star.cutoffs = c(0.05,0.01,0.001))
lrtest(negbin62r, negbin61ar) #significant P-value means model fit
#choose model
lrtest(poisson62r,negbin62r)#significant p-value means we should use negative binomial model

#generate margnal effects
meffectsreturnq <- ggpredict(negbin62r, terms=c("time", "group","product_category")) # generates a tidy data frame  

ggplot(meffectsreturnq,aes(x, predicted, colour=group)) + geom_line(size=1.3) + facet_wrap(~facet)
xlab("bops_in_effect") + ylab("returnquantity") +
  labs(colour="Payment") + 
  scale_colour_discrete(labels=c("Product Category1", "Product Category")) +
  scale_x_continuous(breaks=c(0,1), labels=c("no bops", "bops")) +
  theme(axis.title.x=element_blank())

#check for heterokadasticity
gqtest(negbin62r) 
bptest(negbin62r)
HWrobstder <- sqrt(diag(vcovHC(negbin62r, type="HC1"))) # produces Huber-White robust standard errors 

stargazer(negbin62r,negbin62r,
          se=list(NULL, HWrobstder),
          apply.coef = exp, t.auto=F, p.auto = F,
          title="Regression Results", type="text", 
          column.labels=c("Normal SE", "HW-Robust SE"),
          df=FALSE, digits=3, star.cutoffs = c(0.05,0.01,0.001))  # 


