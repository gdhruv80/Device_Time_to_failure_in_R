# Loading necessary libraries

library(ggplot2)
library(sqldf)
library(survival)
library(survminer)
library(dplyr)
library(OptimalCutpoints)
library(pROC)
library(SDMTools)
library(lift)
library(Amelia)
library(mlbench)
library(corrplot)
library(outliers)

# Overall Initiall Thoughts (Based on Scrolling through the Data and Filtering and Sorting in R)

# 1. The data is time based machine log data for 1168 unique devices for the year 2015

# 2. Would term this as a Reliability Analysis problem very similar to modelling Customer Churn/Customer Adoption
#    time to death/birth of child. Basically its a time to event analysis. Considering the type of problem
#    I chose a more statistical approach to model the survival time of a device and hence used R
#    which has better statistical packages to analyze the same.

# 3. Some initial things that I want to check before I further decide on my modelling technique 
#     (Initiall EDA thoughts) :
     # Is the data censored? If so what type
     # What is the Total failure % for all devices in the dataset. Is it very low (<2-3%) if so 
     # we may need to do resampling.
     # What does the lifecycle of a device look like - Is the hazard rate really high in the begining
     # and end and kind on low during the most of the life in the middle just like say with a bulb/hardisk. If there is any such 
     # pattern evident in the Hazard function we would be better off using paramteric hazard models.
     # Are there gaps between dates when values are recorded for a device and if so we will need to create a new timeline variable
     # to account for that
     # Are there missing values in the dataset

# Loading the dataset
Device_dat <- read.csv("device_failure.csv")

### Basic Descriptive analysis on the dataset -

#Checking Data's overall Time Span

# Formating to required data types

Device_dat["date"] = as.Date(Device_dat[,"date"], format = "%m/%d/%Y")
summary(Device_dat)

# Date Time Span 
# min - 2015-01-01
# max - 2015-11-02

# Checking for missing data

apply(is.na(Device_dat),2,sum)
# No missing data anywhere
# Not required - #missmap(Device_dat[1:300,], col=c("black", "grey"), legend=FALSE)


# Percent of failure cases in the dataset

Device_dat_grp =  aggregate(failure ~ device, Device_dat, max)

fail = (sum(Device_dat_grp["failure"])/nrow(Device_dat_grp))*100
fail
# No of unique devices - 1168
# No of cases of failure - 106 - 9.07% failure rate/ positive cases in the set - Not as low as
# I initially thought so no major concerns here

# Understanding Failure Patterns - Creating a Failure Timeline 

# Creating a easy Timeline variable for each device since its inception - Assuming for now that there
# are no missing days between start and end
Device_dat$Nth_Day_frm_Begin <- ave(as.numeric(Device_dat$date,origin="1970-01-01"), 
                                 Device_dat$device, FUN = seq_along)

# Now Checking if records for a device are missing days in between 
Continuity_chk =  aggregate( date~ device, Device_dat, 
                             FUN = function(x) (max(x)-min(x))/(length(x)-1))

# Counting no of Continous Cases
cat("The % of devices without gap in observation dates is :", (nrow(subset(Continuity_chk, date == 1))/nrow(Continuity_chk))*100)

# As we can see there are roughly 15% of the devices where the timeline is not continous. 
# Hence will have to create a new timeline variable to account for that.


# As 85% is not a very low number therefore out of curosity looking at the timeline for failure 
# rate to see if there is a visible failure pattern.

# Creating the Timeline
s1 = subset(Continuity_chk, date == 1)
s2 = subset(Continuity_chk, date != 1)
Device_Visual <- Device_dat[Device_dat$device %in% s1[,c(1)],]
Device_Visual = Device_Visual[,c('device','Nth_Day_frm_Begin','failure')]

#Note - 1. For the timeline making sure we have significant sample size for each day we dont 
#       want unnecessary spikes in failure percent for a day because of small sample size
#       (excluding Days which have <30 devices recorded on them).
#       2. Also aggregating the data at a weekly level to remove too much noise


Device_Visual_w = Device_dat[c(1:3)]
# Extracting week
Device_Visual_w$week = format(as.Date(Device_Visual_w$date), "%m") 

# Checking if the bulb was tested in that week and if so did it fail or not
# Assumption for creating timeline visual - 
# If the bulb is tested on any day of the week it is assumed to be good for the week

Device_Visual_w_agg =  sqldf("select device,week,
      max(failure) `Failed_or_not`, 
                         count(failure) `No_times_eval_in_week` 
                         from `Device_Visual_w` 
                         group by device,week")

Device_Visual_w_agg_2 =  sqldf("select week,
      sum(Failed_or_not) `Count_Failed`, 
                             count(DISTINCT device) `No_devices_in_week` 
                             from `Device_Visual_w_agg` 
                             group by week")

Device_Visual_w_agg_2$weekly_Failure_Rate = ((Device_Visual_w_agg_2$Count_Failed)/
                                             (Device_Visual_w_agg_2$No_devices_in_week))*100

# Plotting the weekly fail rate using ggplot2
p2 = ggplot(data=Device_Visual_w_agg_2, aes(week, weekly_Failure_Rate)) + 
  geom_col( fill="steelblue") + 
  geom_text(aes(label=round(weekly_Failure_Rate,2)), vjust=1.6, color="white", size=3.5)
 
p2 <- p2 + scale_y_continuous(breaks=seq(0,5,.5))
p2 <- p2 + labs(x = "Weeks after first being Tested",y = "Failure Percent",title = "Failure Timeline For a Device")
p2

# From the weekly hazard rate we cannot see any strong patterns in hazard rates such as 
# a constant/declining/high-constant-high hazard rate. In fact it is pretty erratic. Therefore
# I would prefer to not use a weibull/exp regression or any other parametric model but would rather
# prefer using something like Cox - Regression which does not force any assumptions on the shape of
# underlying data


# Populating correct nth Day from Begining - (Accounting for time gaps in recordings for a device)

# The idea while creating this synthetic variable is to figure out rows where the sucessive date difference
# for a device is greater than 1 and then add that difference to the older date to get the newer one.

# Creating a empty DF
df = Device_dat[FALSE,]
df$Nth_Day_frm_Begin_2 <- df$Nth_Day_frm_Begin

uniq = unique(Device_dat$device)


for (i in 1:length(uniq)){
  # Subsetting the Dataset for each device
  d = subset(Device_dat, device == uniq[i])
  d = d[order(d$date),]
  
  d$Nth_Day_frm_Begin_2[1] = 1
  for (i in 1:(nrow(d)-1)){
    a = d$date[i+1]-d$date[i]
    # For each device correctly populating the timeline to get it from the date level to 
    # a lifecycle (Day1 ... Day 200..) level
    if(a==1){
      d$Nth_Day_frm_Begin_2[i+1] = d$Nth_Day_frm_Begin_2[i] + 1
    }else{  d$Nth_Day_frm_Begin_2[i+1] = a +  d$Nth_Day_frm_Begin_2[i] }
  }
  # Merging all the individual devices together
  df = rbind(df,d)
}

# Another interesting observations here is that the covariates are changing with time. Ie we
# will have to use cox regression with time varying covariates. This would require some more
# data preparation. We will now have to keep all the rows for a device wherever there is a change
# in any of the attributes.

# Removing redundant rows for each device (ie No changes in Attributes)
df_1 = df[FALSE,]
for (i in 1:length(uniq)){
  d = subset(df, device == uniq[i])
  d = d[with(d, c(TRUE, diff(as.numeric(interaction(attribute1, attribute2, attribute3, 
                                                    attribute4, attribute5, attribute6, 
                                                    attribute7, attribute8, attribute9))) != 0)), ]
  df_1 = rbind(df_1,d)
}
df_1[1:10,]
df_2 = df_1

# Further we will have to create 2 new variables T0 and T1 which reflect the begin time and end time for each
# row in the dataset for the device
# Create T0 and T1 - As per documentation from surv package 
# Note that the Time interval for each row is (T0-T1]

df_3 = df_2[FALSE,]
for (i in 1:length(uniq)){
  d = subset(df_2, device == uniq[i])
  d = d[order(d$Nth_Day_frm_Begin_2),]
  
  d$T0[1] = 0 #For each device the begin time for first recording will be 0 
  # Populating T0 and T1 for all rows of a device
  for (i in 1:(nrow(d)-1)){
    d$T1[i] = d$Nth_Day_frm_Begin_2[i+1] - 1
    d$T0[i+1] = d$T1[i]
  }
  
  d$T1[nrow(d)] = d$Nth_Day_frm_Begin_2[nrow(d)]
  # Merging all the devices together
  df_3 = rbind(df_3,d)
}
  
# Removing redundant columns before modelling  
df_4 = subset(df_3, select = -c(Nth_Day_frm_Begin,Nth_Day_frm_Begin_2,date))

df_4[1:10,]
df_xx = subset(df_3,device == 'S1F01R2B')

# Initial Data Analysis - 
# Outlier Detection/Testing - Looking at correlations

# calculate correlations matrix
correlations <- cor(df_4[,3:11])
# create correlation plotpar(mfrow=c(1,5))
corrplot(correlations, method="circle")

# Variable 7 and 8 have a 100% corellation. Infact they are identical. Attribute 3 and 9 have
# significant +ve correlation(good to know). There are no negetive correlations among the variables. 


# Looking at boxplots for all attributes to identify oultiers

# Attribute 1-5 Boxplots
par(mfrow=c(1,5))
for(i in 3:7) {
  boxplot(df_4[,i], main=names(df_4)[i])
}

# Attribute 6-9 Boxplots
par(mfrow=c(1,4))
for(i in 8:11) {
  boxplot(df_4[,i], main=names(df_4)[i])
}

# It is interesting to see that for 6 of the 9 variables more than 75% of the values are 0
# hence rendering whiskers in boxplots kind of useless

# Subset 1 - Variables with boxplots with non zero whiskers
# Checking the % of data outside the whiskers for these 
for(i in c(3,7,8)){
  o = boxplot(df_4[,i], range = 2,plot=FALSE)$out
  print ((length(o)/nrow(df_4))*100)
}


# For these variable Capping outliers to the 1.5 IQR value 

for(i in c(3,7,8)){
  
  # Get the Min/Max values
  max_1 <- quantile(df_4[,i],0.75, na.rm=TRUE) + (IQR(df_4[,i], na.rm=TRUE) * 2 )
  min_1 <- quantile(df_4[,i],0.25, na.rm=TRUE) - (IQR(df_4[,i], na.rm=TRUE) * 2 )
  
  for(s in 1:nrow(df_4)){
    if(df_4[s,i]<min_1){df_4[s,i]=min_1} else if(df_4[s,i]>max_1){df_4[s,i]=max_1} 
  }
}                           

#For attributes 2,3,4,7,9 with boxplots with no whiskers Visual inspection shows that
# all of these have very few in number ~0.001% but significant in magnitude outliers
# Hence Capping extreme outliers to the 0.001%ile value

par(mfrow=c(1,5))
for(i in c(4,5,6,9,11)) {
  boxplot(df_4[,i], main=names(df_4)[i])
}

for(i in c(4,5,6,9,11)){
  
  # Get the Min/Max values
  max_v <- quantile(df_4[,i],0.999, na.rm=TRUE)
  
  for(s in 1:nrow(df_4)){
    if(df_4[s,i]>max_v){df_4[s,i]=max_v} 
  }
}          

# Looking at all the attributes again to see how they look after outlier handling
par(mfrow=c(1,5))
for(i in 3:7) {
  boxplot(df_4[,i], main=names(df_4)[i])
}
par(mfrow=c(1,4))
for(i in 8:11) {
  boxplot(df_4[,i], main=names(df_4)[i])
}

df_5 = df_4

# Removing the attribute 8 which is redundant
df_5 = df_5[,c(-10)]

# MODEL DEVELOPMENT
# Dividing into Train Test and Validate() - Train (65%) Test (20%) Validate (15%)
# Note :- Final declared model results would be the ones obtained from the Validate Set

set.seed(104)
u1 = as.data.frame(sample(uniq))

# Figuring out Unique devices in each test,train,validate
train = sample_frac(u1, 0.65)
#759
test = sample_frac(anti_join(u1, train, by = "sample(uniq)"),0.571)
#234
val = anti_join(u1, rbind(train,test), by = "sample(uniq)")
#175

# Breaking it up in the 3 (Test,Train,Validate)
train_d = subset(df_5, device %in% train$`sample(uniq)`)
test_d = subset(df_5, device %in% test$`sample(uniq)`)
val_d = subset(df_5, device %in% val$`sample(uniq)`)

# Nrow check
nrow(train_d) + nrow(test_d) + nrow(val_d) == nrow(df_4)

# Checking failure cases ratio in all three sets - They are not dispropotionate
a = list(train_d,test_d,val_d)
for (i in a){
  print ((sum(i$failure)/length(unique(i$device)))*100)
}


# Training the Model on the Train Set

# Note - First I have created simple single variable COX models for each of the attribute
# and selected only those variables to build the final model which have significant/non zero 
# coffecients (which is checked by looking at the Pvalue for each of these variables)
# As there are no significant correlations between any of the input variables at this stage
# this will be a robust technique


# Spitting out significant variables

for (i in c(3:10)){
  if(i==3){print ("Independent Variables with significant coffecients P<0.05 :-")}
res.cox_i =  coxph(Surv(T0,T1, failure) ~train_d[,i] , data = train_d)
res.cox_i_2 = data.frame(summary(res.cox_i)$coefficients)
if(res.cox_i_2[1,5]<0.05){cat ("attribute",(i-2),"\n")}
}

# Building a final model based on the 3 attributes selected above
res.cox <- coxph(Surv(T0,T1, failure) ~ attribute2+ attribute4 + attribute7, data = train_d) 
summary(res.cox)
# The model has a R(square) of 0.002 and is significant based on all 3 tests (Log-Liklihood,Wald and Logrank)
# Hence we have a statistically valid model with significant coffecients that we can go forward with

# Looking at the survival plot to see the overall survivial probability (assuming mean values for all Covariates)
# and CI's to get a feel of the model results

ggsurvplot(survfit(res.cox), color = "#2E9FDF",
           ggtheme = theme_minimal())


# Understanding the importance/effects of covariates on the survival probabilities

#Example - 1
# Seeing how doubling up the value of attribute 2 affects the survival probability
HR_attribute_2 <- with(train_d,data.frame(attribute2 = c(10000, 20000), 
                          attribute4 = mean(attribute4, na.rm = TRUE),
                          attribute7 = mean(attribute7, na.rm = TRUE)))

res.cox_hr = survfit(res.cox, newdata = HR_attribute_2)
ggsurvplot(res.cox_hr, conf.int = TRUE, legend.labs=c("Attrib_2=10000", "Attrib_2=20000"),
           ggtheme = theme_minimal())

# Very interesting to see that on doubling the covariate value for attribute 2 the probability of survivial 
# at the 300 day mark goes down from 59% to 40%. Hence describing the positive hazard ratio
# that we obtained for this variable before.
# Similar plots can be made for variable 4 and 7

# Running DAIGNOSTICS ON COX REGRESSION

# It is important to check if the assumptions for the cox model are met by the data or not

# 1.
# Checking propotional hazard assumption (if HR is constant over time or not)
Prop_check_1 = cox.zph(res.cox)
Prop_check_1


# rho chisq      p
# attribute2 0.253  4.45 0.0350
# attribute4 0.156  3.93 0.0474
# attribute7 0.112  1.05 0.3059
# GLOBAL        NA 10.68 0.0136


# Using Schenfeld residuals to visually double check the p values obtained above
Prop_check_graph = ggcoxzph(Prop_check_1)
Prop_check_graph

# Note - It is interesting to see that for both attribute 2 and 4 the assumptions that hazards
# are propotional over time or hazard ratios are constant are not valid. P>0.05 indicating that
# residuals do show a strong pattern over time. This can also be verified by looking at the Schenfeld 
# residual plots. Hence we will have to recreate th model by including interaction terms of 
# attribute 2 and 4 with time (T0). This should accomodate for the changing hazard ratios.

# Cox Model version 2 - Accounting for non proportional hazards 

res.cox_2 = coxph(Surv(T0,T1, failure) ~ 
                    attribute2:T0+ attribute4:T0 + attribute7, data = train_d) 

summary(res.cox_2)
# Similar R(square) with decreases in the hazarad ratios for interaction attributes 2 and 4 as expected.
# and a statistically significant model

# 2.
#Indentifying if any of data points significantly affect the regression coffeccient values.

# Looking at the max change in Beta values for influential data points.

ggcoxdiagnostics(res.cox_2, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
# As expected as we already accounted for extreme outliers before none of the data points individually
# affects the Beta values in a significant/drastic fashion


#After completing Diagnostics & Training Moving Forward with making predictions on the Test Set

# Predicting on test set and Computing ROC curve and AUC

test_pred_c = predict(res.cox_2, test_d,
        type=c("lp"),collapse =test_d$device)

score_test = data.frame(test_pred_c)

# Merging predictions with X values
score_test_1 = merge(score_test,
                     sqldf("select device,max(failure) `Failed_or_not`, count(failure) 'days of data'
                           from `test_d` group by device"),by.x = 0,by.y = "device")      

              
# Plotting ROC curve 
plot(roc(score_test_1$Failed_or_not, score_test_1$test_pred_c, direction="<"),
     col="yellow", lwd=3, main="The turtle finds its way")              

# Area under ROc = 0.723 - Considering we just had 8 independent covariates thats not a bad AUC

# Minimizing False Positives while Maximizing True Positives at the same time

# Getting optimal operation/threshold point (using Youden criteria: Sensitivity + Specificity - 1)
# Now as the costs for False Positive and False Negetives were not provided and depend on the business scenario
# I have assumed them both to be 1. But the cutpoint criteria below provides us with the option of obtaining
# a new cutoff point based on different costs for these.
              
optimal.cutpoint.Youden <-optimal.cutpoints(costs.ratio = 1,CFP = 1, CFN = 1,X = "test_pred_c", 
                                            status = "Failed_or_not",tag.healthy = 0,
                                            methods = "Youden", data = score_test_1)

# Getting the confusion matrix based on our model operation cutoff

CF = confusion.matrix(score_test_1$Failed_or_not,score_test_1$test_pred_c,
                      threshold=optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$cutoff)

CF
cat ("For the Test Set The True Positive Rate is :",
    CF[4]/sum(CF[3],CF[4]), "while the False Positive Rate is:",round(CF[2]/sum(CF[2],CF[1]),2))


## Verfying Test set predictions by Predicting on validation set 


val_pred_c = predict(res.cox, val_d,
                      type=c("lp"),collapse =val_d$device);

score_val = data.frame(val_pred_c);

score_val_1 = merge(score_val,
                     sqldf("select device,max(failure) `Failed_or_not`, count(failure) 'days of data'
                           from `val_d` group by device"),by.x = 0,by.y = "device")   

# ROC
plot(roc(score_val_1$Failed_or_not, score_val_1$val_pred_c, direction="<"),
     col="yellow", lwd=3, main="The turtle finds its way")                     

# Validation ROC looks pretty similar to test set indicating that we were not overfitting 
# on the test set

# Using the cutoff from the Test set to get the confusion matrix        

CF_V = confusion.matrix(score_val_1$Failed_or_not,score_val_1$val_pred_c,
                        threshold=optimal.cutpoint.Youden$Youden$Global$optimal.cutoff$cutoff)

CF_V             
cat ("For the Test Set The True Positive Rate is :",
     round(CF_V[4]/sum(CF_V[3],CF_V[4]),2), "while the False Positive Rate is:",round(CF_V[2]/sum(CF_V[2],CF_V[1]),2))

# very similar results here as the Test set. So finally with our model we can expect to
# correctly identify ~60% of the failures before they happen while just raising a false alarm
# for 8% of the cases.

#Lift calculation for the validation set
#Note : Lift calculation is a very important calculation here especially because the total
# number of positive (failure) cases in our model are so small and hence just providing the AUC could be misleading. 
# It really enables us to show the strength of our model. The overall lift of 4.7 here shows that our model has a 4.7X
# fold/time capability of identifying failures before hapenning(59%) as compared to if there was no
# model in place and we were maintaing all devices (just 12.5%) while just raising false errors for only 8% 
# of the cases that it classifies.

Errors_pnct = sum(CF_V[3],CF_V[4])/sum(CF_V)
Lift = (CF_V[4]/sum(CF_V[3],CF_V[4]))/Errors_pnct
cat ("The lift provided by the Model is:",Lift)


# Next Steps :- 

# As next steps I would like to definitly apply some more models here and compare 
# performace. Would Even experiment with some parametric models and see how they scores against the
# cox. Also what would be interesting is to see if gathering more data esentially specific data
# would help us improve our model performance.

# 1 simple test to see if collecting more data especially for devices which are in our dataset
# for shorter span(right censored) would help.

# Trying to see if the average days that we have in the dataset differ for devices those which the 
# model predicted correctly vs those it misclassified. 
score_val_1$Pred_1_0 = sapply(score_val_1$val_pred_c,
                              function(x){if(x>0.58){y =1}else{y = 0}})

score_val_1$Miss_class =  ifelse(score_val_1$Pred_1_0 == score_val_1$Failed_or_not,0,1)

print (aggregate(`days of data` ~ Miss_class, score_val_1, mean))

# The average no of days are pretty close a difference could have indicated that we are 
# misclassifing those records for which we have less data/no of days of info hence it 
# would then make sense to collect more data for these devices.

