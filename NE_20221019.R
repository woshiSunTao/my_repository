rm(list = ls())

library(car)
library(ggplot2)
library(glmnet)
library(MASS)
library(survival)
library(rms)
library(survminer)
library(ggridges)
library(pROC)
library(plotROC)
library(riskRegression)
library(magrittr)
library(DynNom)
library(packrat)
library(rsconnect)
library(readxl)
library(plyr)
library(survivalROC)
library("gtsummary")
library(dplyr)
library(flextable)
library(Hmisc)
library(mlr)
library(magicfor)
library(readxl)
library(table1)
library(survival)
library(survivalROC)
library(survminer)
library(rms)
library(magicfor)
library(randomForestSRC)
library(pec)
library(mlr)
library(partykit)
library(party)

data1 <- read.csv("20220919.csv")

train <- data1
#train=subset(train,train$Histologic.Type.ICD.O.3==8041|train$Histologic.Type.ICD.O.3==8045)

train=subset(train,train$Histologic.Type.ICD.O.3==8260|train$Histologic.Type.ICD.O.3==8050)

###SCC: N=188624 ICD:8070,8071,8072,8083
train=subset(train,train$Diagnostic.Confirmation=="Positive histology")
###确诊方式为病理活检：n=24830  N=163794
#  train=subset(train,train$Sequence.number=="One primary only")
###One primary only，n=53694,N=110100

train <- subset(train,train$Year.of.diagnosis>2003)

### T stage Unknown n=163 N=68010
table(train$Primary.Site...labeled)
#train=subset(train,train$Primary.Site...labeled!="C34.9-Lung, NOS")
train$Primary_site=ifelse(train$Primary.Site...labeled=="C34.0-Main bronchus"|train$Primary.Site...labeled=="C34.8-Overlapping lesion of lung","Other",train$Primary.Site...labeled)
train$Primary_site=ifelse(train$Primary.Site...labeled=="C34.1-Upper lobe, lung","Upper",train$Primary_site)
train$Primary_site=ifelse(train$Primary.Site...labeled=="C34.2-Middle lobe, lung","Middle",train$Primary_site)
train$Primary_site=ifelse(train$Primary.Site...labeled=="C34.3-Lower lobe, lung","Lower",train$Primary_site)
train$Primary_site=ifelse(train$Primary.Site...labeled=="C34.9-Lung, NOS","L_NOS",train$Primary_site)


table(train$Year.of.diagnosis)
nrow(train)

train$Tumor_size=ifelse(train$CS.tumor.size..2004.2015.=="Blank(s)",train$Tumor.Size.Summary..2016..,train$CS.tumor.size..2004.2015.)
train=subset(train,train$Tumor_size!="Blank(s)")
train$Tumor_size=as.numeric(train$Tumor_size)

train=subset(train,train$Tumor_size!=0)

train=subset(train,train$Tumor_size!=990)
train=subset(train,train$Tumor_size!=991)
train=subset(train,train$Tumor_size!=992)
train=subset(train,train$Tumor_size!=993)
train=subset(train,train$Tumor_size!=994)
train=subset(train,train$Tumor_size!=995)
train=subset(train,train$Tumor_size!=996)
train=subset(train,train$Tumor_size!=997)
train=subset(train,train$Tumor_size!=998)


nrow(train)


train$N_stage=train$Derived.AJCC.N..6th.ed..2004.2015.
train$N_stage=ifelse(train$N_stage=="Blank(s)"|train$N_stage=="NX",train$Derived.SEER.Combined.N..2016.2017.,train$N_stage)
train$N_stage=ifelse(train$N_stage=="Blank(s)"|train$N_stage=="NX",train$Derived.EOD.2018.N..2018..,train$N_stage)
train=subset(train,train$N_stage!="Blank(s)"&train$N_stage!="88"&train$N_stage!="NX"&train$N_stage!="pX"&train$N_stage!="cX")
train$N_stage=ifelse(train$N_stage=="c0","N0",
                     ifelse(train$N_stage=="c1","N1",
                     ifelse(train$N_stage=="c2","N2",
                     ifelse(train$N_stage=="c3","N3",
                     ifelse(train$N_stag=="p0","N0",
                     ifelse(train$N_stage=="p1","N1",
                     ifelse(train$N_stage=="p2","N2",
                     ifelse(train$N_stage=="p3","N3",train$N_stage))))))))

table(train$CS.extension..2004.2015.)
train$adjustedSixthMstage=ifelse(train$CS.extension..2004.2015.=="720","M1",train$Derived.AJCC.M..6th.ed..2004.2015.)
train$adjustedSixthMstage=ifelse(train$CS.extension..2004.2015.=="760","M1",train$Derived.AJCC.M..6th.ed..2004.2015.)
train$adjustedSixthMstage=ifelse(train$CS.extension..2004.2015.=="790","M1",train$Derived.AJCC.M..6th.ed..2004.2015.)
train$M_stage=ifelse(train$Derived.EOD.2018.M..2018..=="Blank(s)",train$Derived.SEER.Combined.M..2016.2017.,train$Derived.EOD.2018.M..2018..)
train$M_stage=ifelse(train$M_stage=="Blank(s)",train$Derived.AJCC.M..7th.ed..2010.2015.,train$M_stage)
train$M_stage=ifelse(train$M_stage=="Blank(s)",train$adjustedSixthMstage,train$M_stage)
train=subset(train,train$M_stage!="M1NOS"&train$M_stage!="MX")
library(tidyverse)
train$M_stage=ifelse(str_detect(train$M_stage,"0"),"M0","M1")
table(train$M_stage)
###M stage Unknown n=376  N=67334


###N stage Unknown n=2091 N=71026
#train=subset(train,train$Combined.Summary.Stage..2004..!="Unknown/unstaged")
#train$M_stage=ifelse(train$Combined.Summary.Stage..2004..=="Distant","Yes","No")
###M stage Unknown n=208 N=70818
train$T_stageAdjustedSixth=train$Derived.AJCC.T..6th.ed..2004.2015.
table(train$T_stageAdjustedSixth,train$Year.of.diagnosis)
train$T_stageAdjustedSixth=ifelse(train$T_stageAdjustedSixth=="550"&train$T_stageAdjustedSixth=="T3","T2",train$T_stageAdjustedSixth)
train$T_stageAdjustedSixth=ifelse(train$T_stageAdjustedSixth=="500"&train$T_stageAdjustedSixth=="T3","T2",train$T_stageAdjustedSixth)
train$T_stageAdjustedSixth=ifelse(train$CS.extension..2004.2015.=="650"&train$T_stageAdjustedSixth=="T4","T3",train$T_stageAdjustedSixth)
train$T_stageAdjustedSixth=ifelse(train$T_stageAdjustedSixth=="T2"&train$Tumor_size>50&train$Tumor_size<71,"T3",train$T_stageAdjustedSixth)
train$T_stageAdjustedSixth=ifelse(train$CS.site.specific.factor.1..2004.2017.varying.by.schema.=="10","T3",train$T_stageAdjustedSixth)
train$T_stageAdjustedSixth=ifelse(train$CS.site.specific.factor.1..2004.2017.varying.by.schema.=="10","T3",train$T_stageAdjustedSixth)
train$T_stageAdjustedSixth=ifelse(train$T_stageAdjustedSixth=="T4","T4",train$T_stageAdjustedSixth)
train$T_stageAdjustedSixth=ifelse(train$T_stageAdjustedSixth=="T3"&train$Tumor_size>70,"T4",train$T_stageAdjustedSixth)
train$T_stageAdjustedSixth=ifelse(train$T_stageAdjustedSixth=="T2"&train$Tumor_size>70,"T4",train$T_stageAdjustedSixth)
train$T_stageAdjustedSeventh=ifelse(train$Derived.AJCC.T..7th.ed..2010.2015.=="Blank(s)",train$Derived.SEER.Combined.T..2016.2017.,train$Derived.AJCC.T..7th.ed..2010.2015.)
train$T_stageAdjustedSeventh=ifelse(train$T_stageAdjustedSeventh=="T1a"|train$T_stageAdjustedSeventh=="T1b"|train$T_stageAdjustedSeventh=="T1NOS","T1",train$T_stageAdjustedSeventh)
train$T_stageAdjustedSeventh=ifelse(train$T_stageAdjustedSeventh=="T2a"|train$T_stageAdjustedSeventh=="T2b"|train$T_stageAdjustedSeventh=="T2NOS","T2",train$T_stageAdjustedSeventh)
train$T_stageAdjustedSeventh=ifelse(train$CS.extension..2004.2015.=="550"&train$T_stageAdjustedSeventh=="T3","T2",train$T_stageAdjustedSeventh)
train$T_stageAdjustedSeventh=ifelse(train$CS.extension..2004.2015.=="500"&train$T_stageAdjustedSeventh=="T3","T2",train$T_stageAdjustedSeventh)
train$T_stageAdjustedSeventh=ifelse(train$Tumor_size>50&train$Tumor_size<71&train$T_stageAdjustedSeventh=="T2","T3",train$T_stageAdjustedSeventh)
train$T_stageAdjustedSeventh=ifelse(train$Tumor_size>70&train$T_stageAdjustedSeventh=="T3","T4",train$T_stageAdjustedSeventh)
train$Derived.EOD.2018.T..2018..=ifelse(train$Derived.EOD.2018.T..2018..=="T1a"|train$Derived.EOD.2018.T..2018..=="T1b"|train$Derived.EOD.2018.T..2018..=="T1c"|train$Derived.EOD.2018.T..2018..=="T1mi","T1",train$Derived.EOD.2018.T..2018..)
train$Derived.EOD.2018.T..2018..=ifelse(train$Derived.EOD.2018.T..2018..=="T2a"|train$Derived.EOD.2018.T..2018..=="T2b","T2",train$Derived.EOD.2018.T..2018..)
train$T_stage=ifelse(train$Derived.EOD.2018.T..2018..=="Blank(s)"|train$Derived.EOD.2018.T..2018..=="TX",train$T_stageAdjustedSeventh,train$Derived.EOD.2018.T..2018..)
train$T_stage=ifelse(train$T_stage=="Blank(s)"|train$T_stage=="TX",train$T_stageAdjustedSixth,train$T_stage)
nrow(train)
library(stringr)
train$T_stage=ifelse(str_detect(train$T_stage,"1"),"T1",train$T_stage)
train$T_stage=ifelse(str_detect(train$T_stage,"2"),"T2",train$T_stage)
train$T_stage=ifelse(str_detect(train$T_stage,"3"),"T3",train$T_stage)
train$T_stage=ifelse(str_detect(train$T_stage,"4"),"T4",train$T_stage)
train=subset(train,train$T_stage!="c0"&train$T_stage!="cX"&train$T_stage!="Tis"&train$T_stage!="T0"&
               train$T_stage!="Blank(s)"&train$T_stage!="TX")
### T stage Unknown n=1023 N=70004
nrow(train)

#train <- subset(train,train$M_stage=="M0")
#train <- subset(train,train$N_stage=="N0"|train$N_stage=="N1")
#train <- subset(train,train$T_stage=="T1"|T_stage=="T2"|T_stage=="T3")

#library(tidyr)
#train <- tidyr::unite(train,"NT",T_stage,N_stage,sep="",remove=FALSE)

#train <- subset(train,train$NT != "T3N1")


### Primary site Unknown n=2645  N=68010
train=subset(train,train$Laterality=="Left - origin of primary"|train$Laterality=="Right - origin of primary")
###Laterality Unknown n=179  N=67831

library(tidyverse)
train$Age <- train$Age.recode.with.single.ages.and.100.
train <- subset(train,train$Age != "100+ years")
train <- separate(data = train, col = Age, into = c("Age", "a1"), sep = " ")
train$Age <- as.numeric(train$Age) 
train <- subset(train,train$Age.recode.with.single.ages.and.100. >= 18)

train <- subset(train,train$Race.recode..White..Black..Other. != "Unknown")

train <- subset(train,(train$Survival.months) != "Unknown")


train$Race <- train$Race.recode..White..Black..Other.
train$Surgery <- train$RX.Summ..Surg.Prim.Site..1998..
train$Radiation <- train$Radiation.recode
train$Chemotherapy <- train$Chemotherapy.recode..yes..no.unk.
train$Histologic <- train$Histologic.Type.ICD.O.3

train$Race <- ifelse(train$Race=="Black","Black",ifelse(train$Race=="White","White","Other"))
train$Laterality <- ifelse(train$Laterality=="Left - origin of primary","Left","Right")
nrow(train)

train <- subset(train,train$Surgery != 99)
train$Surgery <- ifelse(train$Surgery==0,"No_Surgery",ifelse(train$Surgery>=15&train$Surgery<=25,"Sublobectomy",ifelse(train$Surgery>=30&train$Surgery<=50,"Lobectomy",ifelse(train>=55&train<=70,"Pnermonectomy","Other_Surgery"))))

nrow(train)

train$Marital <- ifelse(train$Marital.status.at.diagnosis=="Married (including common law)","Married","Other")

train$Radiation <- ifelse(train$Radiation.recode=="None/Unknown"|train$Radiation.recode=="Recommended, unknown if administered"|train$Radiation.recode=="Refused (1988+)","No/unknown","Yes")

train$Grade <- ifelse(train$Grade..thru.2017.=="Moderately differentiated; Grade II"|train$Grade..thru.2017.=="Well differentiated; Grade I","I-II",ifelse(train$Grade..thru.2017.=="Poorly differentiated; Grade III"|train$Grade..thru.2017.=="Undifferentiated; anaplastic; Grade IV","III-IV","Unknown"))

train$Survivalmonths <- as.numeric(train$Survival.months)
train$status <- ifelse(train$Vital.status.recode..study.cutoff.used. =="Alive",0,1)
train$specific <- ifelse(train$SEER.cause.specific.death.classification=="Dead (attributable to this cancer dx)",1,0)


nrow(train)

raw <- train[,c("Age","Sex","Race","Laterality","T_stage","N_stage","M_stage","Marital","Histologic","Grade","Tumor_size","Primary_site","Chemotherapy","Surgery","Radiation","Survivalmonths","specific","status")]
nrow(raw)

raw$Age <- ifelse(raw$Age<79,"<79",">=79")

raw$Tumor_size <- ifelse(raw$Tumor_size<28,"<28",ifelse(raw$Tumor_size>=28&raw$Tumor_size<=52,"28-52",ifelse(raw$Tumor_size>52&raw$Tumor_size<990,">52","Unknown")))

ordered.item  <- c("<28","28-52",">52","Unknown")
ordered.item <- factor(1:length(ordered.item),labels = ordered.item)
raw$Tumor_size <- factor(raw$Tumor_size,levels = levels(ordered.item))

raw$Surgery_group <- ifelse(raw$Surgery=="Lobectomy","Lobectomy",ifelse(raw$Surgery=="No_Surgery","No_Surgery","Other_Surgery"))
ordered.item  <- c("No_Surgery","Other_Surgery","Lobectomy")
ordered.item <- factor(1:length(ordered.item),labels = ordered.item)
raw$Surgery_group <- factor(raw$Surgery_group,levels = levels(ordered.item))


raw <- raw[,c("Age","Sex","Race","Laterality","T_stage","N_stage","M_stage","Marital","Grade","Tumor_size","Primary_site","Chemotherapy","Surgery_group","Radiation","Survivalmonths","specific")]
nrow(raw)


raw$Age <- factor(raw$Age)
raw$Sex <- factor(raw$Sex)
raw$Race <- factor(raw$Race)
raw$Laterality <- factor(raw$Laterality)
raw$T_stage <- factor(raw$T_stage)
raw$N_stage <- factor(raw$N_stage)
raw$M_stage<- factor(raw$M_stage)
raw$Marital<- factor(raw$Marital)
raw$Grade<- factor(raw$Grade)
raw$Primary_site<- factor(raw$Primary_site)
raw$Chemotherapy <- factor(raw$Chemotherapy)
raw$Radiation <- factor(raw$Radiation)

#write.csv(raw,"20221106_1.csv")


raw$three_y_specific <- ifelse(raw$Survivalmonths<36&raw$specific==1,1,0)

raw$five_y_specific <- ifelse(raw$Survivalmonths<60&raw$specific==1,1,0)

raw$ten_y_specific <- ifelse(raw$Survivalmonths<120&raw$specific==1,1,0)

raw <- subset(raw,raw$Survivalmonths>0)

#############################################################################################
fit <- survfit(Surv(Survivalmonths,specific) ~ 1 , data = raw)
ggsurvplot(fit, 
           data = raw,
           pval = TRUE,
           pval.coord = c(0, 0.03), 
           surv.median.line = "hv", 
           legend.title="Type",
           xlim=c(1,181),
           palette=c("#005824", "#E41A1C","#377EB8","#984EA3"),
           #title="Overall Survival",
           ylab="Cumulative survival (percentage)",xlab = " Survival time (Months)", #???ĺ???????
           censor.shape = 124,censor.size = 2,conf.int = FALSE,
           break.x.by = 30,
           risk.table = TRUE, # 添加风险表
           risk.table.col = "strata", # 根据分层更改风险表颜色
           # xlab = "PDTC CSS(m)", # 指定x轴标签
           ggtheme = theme_bw())







#write.csv(raw,"20221111.csv")


threedata <- raw[,c("Age","Sex","Race","Laterality","T_stage","N_stage","M_stage","Marital","Grade","Tumor_size","Primary_site","Chemotherapy","Surgery_group","Radiation","Survivalmonths","three_y_specific")]

threedata <- subset(threedata,threedata$Survivalmonths != 0)

#write.csv(threedata,"20221106.csv")

#cox model in complete cases
random_tune <- makeTuneControlRandom(maxit = 1L)
rdesc = makeResampleDesc( "CV", iters = 10, stratify = TRUE ) 
###1
seeds<-c(1:50)
magic_for(print, silent = TRUE)
for (k in seeds){
  set.seed(k)
  index <- sample(1:nrow(threedata), round(0.7*nrow(threedata)))
  test_set <- threedata[-index,]
  train_set<- threedata[index,]
  task<-makeSurvTask(data = train_set,target=c('Survivalmonths','three_y_specific'))
  cox.lrn <- makeLearner(cl="surv.coxph", 
                         predict.type="response")
  modcox = train(cox.lrn, task)
  train_pred<-predict(modcox, newdata = train_set)
  train_p<-performance(train_pred, measures = list(cindex)) 
  test_pred<-predict(modcox, newdata = test_set)
  test_p<-performance(test_pred, measures = list(cindex)) 
  print(round(train_p,3),round(test_p,3))
}
performance<-magic_result_as_dataframe()
summary(performance)

#%%%%%%%%%%%%%%%%%%%%% Gradient boosting with smooth components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
random_tune <- makeTuneControlRandom(maxit = 1L)
rdesc = makeResampleDesc( "CV", iters = 10, stratify = TRUE ) 
###1
seeds<-c(1:100)
magic_for(print, silent = TRUE)
for (k in seeds){
  set.seed(k)
  index <- sample(1:nrow(threedata), round(0.3*nrow(threedata)))
  test_set <- threedata[-index,]
  train_set<- threedata[index,]
  task<-makeSurvTask(data = train_set,target=c('Survivalmonths','three_y_specific'))
  cox.lrn <- makeLearner(cl="surv.gamboost", 
                         predict.type="response")
  modcox = train(cox.lrn, task)
  train_pred<-predict(modcox, newdata = train_set)
  train_p<-performance(train_pred, measures = list(cindex)) 
  test_pred<-predict(modcox, newdata = test_set)
  test_p<-performance(test_pred, measures = list(cindex)) 
  train_p;test_p
  
  print(round(train_p,3),round(test_p,3))
}
performance<-magic_result_as_dataframe()
summary(performance)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%% Gradient Boosting Machine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
random_tune <- makeTuneControlRandom(maxit = 1L)
rdesc = makeResampleDesc( "CV", iters = 10, stratify = TRUE ) 
###1
seeds<-c(1:100)
magic_for(print, silent = TRUE)
for (k in seeds){
  set.seed(k)
  index <- sample(1:nrow(threedata), round(0.3*nrow(threedata)))
  test_set <- threedata[-index,]
  train_set<- threedata[index,]
  task<-makeSurvTask(data = train_set,target=c('Survivalmonths','three_y_specific'))
  cox.lrn <- makeLearner(cl="surv.gbm", 
                         predict.type="response")
  modcox = train(cox.lrn, task)
  train_pred<-predict(modcox, newdata = train_set)
  train_p<-performance(train_pred, measures = list(cindex)) 
  test_pred<-predict(modcox, newdata = test_set)
  test_p<-performance(test_pred, measures = list(cindex)) 
  train_p;test_p
  
  print(round(train_p,3),round(test_p,3))
}
performance<-magic_result_as_dataframe()
summary(performance)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%% Gradient Boosting with Componentwise Linear Models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
random_tune <- makeTuneControlRandom(maxit = 1L)
rdesc = makeResampleDesc( "CV", iters = 10, stratify = TRUE ) 

seeds<-c(1:100)
magic_for(print, silent = TRUE)
for (k in seeds){
  set.seed(k)
  index <- sample(1:nrow(threedata), round(0.3*nrow(threedata)))
  test_set <- threedata[-index,]
  train_set<- threedata[index,]
  task<-makeSurvTask(data = train_set,target=c('Survivalmonths','three_y_specific'))
  cox.lrn <- makeLearner(cl="surv.glmboost", 
                         predict.type="response")
  modcox = train(cox.lrn, task)
  train_pred<-predict(modcox, newdata = train_set)
  train_p<-performance(train_pred, measures = list(cindex)) 
  test_pred<-predict(modcox, newdata = test_set)
  test_p<-performance(test_pred, measures = list(cindex)) 
  train_p;test_p
  
  print(round(train_p,3),round(test_p,3))
}
performance<-magic_result_as_dataframe()
summary(performance)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











# surv.rpart in complete cases
magic_for(print, silent = TRUE)
table(threedata$Survivalmonths)
threedata <- subset(threedata,threedata$Survivalmonths>0)
for (k in seeds){
  set.seed(k)
  index <- sample(1:nrow(threedata), round(0.7*nrow(threedata)))
  test_set <- threedata[-index,]
  train_set <- threedata[index,]
  task<-makeSurvTask(data = train_set,target=c('Survivalmonths','three_y_specific'))
  rpart.lrn<-makeLearner(cl='surv.rpart', predict.type = 'response')
  getParamSet("surv.rpart")
  model_params_2<-makeParamSet(
    makeIntegerParam('minsplit', lower=1,upper=20),
    makeIntegerParam('maxdepth',lower=1,upper=30))
  # Tune model to find best performing parameter settings using random search algorithm
  tuned_model_2 <- tuneParams(learner = rpart.lrn,
                              task = task,
                              resampling = rdesc,
                              measures = cindex,       #  performance measure, this can be changed to one or many
                              par.set = model_params_2,
                              control = random_tune,
                              show.info = FALSE)
  # Apply optimal parameters to model
  rpart.lrn <- setHyperPars(learner = rpart.lrn,
                            par.vals = tuned_model_2$x)
  modrpart = train(rpart.lrn, task)
  train_pred_rpart<-predict(modrpart, newdata = train_set)
  train_p<-performance(train_pred_rpart, measures = list(cindex)) # c-index  in training set
  test_pred_rpart<-predict(modrpart, newdata = test_set) #prediction in test set
  test_p<-performance(test_pred_rpart, measures = list(cindex)) # 
  print(round(train_p,3),round(test_p,3))
}
performance<-magic_result_as_dataframe()
summary(performance)


#surv.randomForestSRC in complete cases
magic_for(print, silent = TRUE)
table(threedata$Sex)
threedata$Sex=ifelse(threedata$Sex=='Female',0,1)

for (k in seeds) {
  set.seed(k)
  index <- sample(1:nrow(threedata), round(0.7*nrow(threedata)))
  train_set <- threedata[index,]
  test_set <- threedata[-index,]
  train_set<-train_set[,c(1:16)]
  test_set<-test_set[,c(1:16)]
  common <- intersect(names(train_set), names(test_set)) 
  for (p in common) { 
    if (class(train_set[[p]]) == "factor") { 
      levels(test_set[[p]]) <- levels(train_set[[p]]) 
      levels(train_set[[p]]) <- levels(test_set[[p]]) 
      # print(levels(test_set[[p]]))
    } 
  }
  task<-makeSurvTask(data = train_set,
                     target=c('Survivalmonths','three_y_specific'))
  rfsrc.lrn<-makeLearner(cl='surv.randomForestSRC',
                         predict.type = 'response')
  getParamSet("surv.randomForestSRC")
  model_params_3<-makeParamSet(
    makeIntegerParam('ntree', lower=10,upper=100),
    makeIntegerParam('mtry',lower = 1,upper = 3),
    #makeIntegerParam('nodedepth',lower = 1, upper = 20),
    #makeIntegerParam('nodesize',lower = 3, upper=100),
    makeIntegerParam('nsplit',lower = 0, upper=10),
    makeDiscreteParam('splitrule',values = 'logrank', special.vals = list('logrank','logrankscore','random'))
  )
  # Tune model to find best performing parameter settings using random search algorithm
  tuned_model_3 <- tuneParams(learner = rfsrc.lrn,
                              task = task,
                              resampling = rdesc,
                              measures =  cindex,       #  performance measure, this can be changed to one or many
                              par.set = model_params_3,
                              control = random_tune,
                              show.info = FALSE)
  # Apply optimal parameters to model
  rfsrc.lrn <- setHyperPars(learner = rfsrc.lrn,
                            par.vals = tuned_model_3$x)
  modrfsrc = train(rfsrc.lrn, task)
  train_pred_rfsrc<-predict(modrfsrc, newdata = train_set)
  train_p<-performance(train_pred_rfsrc, measures = list(cindex)) # c-index  in training set
  test_pred_rfsrc<-predict(modrfsrc, newdata = test_set) #prediction in test set
  test_p<-performance(test_pred_rfsrc, measures = list(cindex)) #  
  print(round(train_p,3),round(test_p,3))
}
performance<-magic_result_as_dataframe()
summary(performance)



lrn = makeLearner("classif.rpart", predict.type = "prob")
mod = train(lrn, task = sonar.task)
pred = predict(mod, task = sonar.task)
cal = generateCalibrationData(pred)
cal$proportion

cal = generateCalibrationData(pred, groups = 3)
cal$proportion
##      Learner           bin Class Proportion
## 1 prediction [0.000,0.267)     M 0.08860759
## 2 prediction [0.267,0.925)     M 0.51282051
## 3 prediction [0.925,1.000]     M 0.93333333
plotCalibration(cal)


cal = generateCalibrationData(pred)
plotCalibration(cal, smooth = TRUE)
## `geom_smooth()` using formula 'y ~ x'


