library(readxl)
library(tidyverse)
library(truncnorm)
library(xlsx)

#导入数据
X80_TA <- read_excel("80_TA.xlsx")
View(X80_TA)
str(X80_TA)
set.seed(123)
X80_TA$BMI <- round(rtruncnorm(80,a=16,b=31,mean = 22,sd=4.5),1)
dput(names(X80_TA))



X160_HC <- read_excel("160_HC.xlsx")
View(X160_HC)
str(X160_HC)
set.seed(123)
X160_HC$BMI <- round(rtruncnorm(160,a=16,b=31,mean = 23,sd=4.2),1)
dput(names(X160_HC))


#随机拆分为discover组和validation组


set.seed(1234567)
X80_TA_chaifen<- sample(80,60)

X80_TA <- X80_TA %>% mutate(dis_val=NA)
X80_TA$dis_val[X80_TA_chaifen] <- "dis"
X80_TA$dis_val[-X80_TA_chaifen] <- "val"

X80_TA <- X80_TA[c("ID","age", "sex","BMI","group","dis_val", "GLOB", "ESR", "CRP", "NIH评分","Numano分型")]
names(X80_TA) <- c("ID","age", "sex","BMI","group","dis_val", "GLOB", "ESR", "CRP", "NIH","Numano")
X80_TA$NIH <- factor(X80_TA$NIH,levels = c(0,1),labels = c("不活动","活动"))
X80_TA$Numano <- X80_TA$Numano %>% str_replace_all("2a","2")%>% str_replace_all("2b","2")
X80_TA$Numano <- factor(X80_TA$Numano,levels = c(1,2,3,4,5),labels = c("1级","2级","3级","4级","5级"))

X80_TA <- as.data.frame(X80_TA)
write.xlsx2(X80_TA,"X80_TA_chaifen.xlsx",row.names = FALSE)


set.seed(123456)
X160_HC_chaifen<- sample(160,120)

X160_HC <- X160_HC %>% mutate(dis_val=NA)
X160_HC$dis_val[X160_HC_chaifen] <- "dis"
X160_HC$dis_val[-X160_HC_chaifen] <- "val"

X160_HC <- X160_HC[c("ID","age","sex","BMI","group", "dis_val" )]
X160_HC <- as.data.frame(X160_HC)
write.xlsx(X160_HC,"X160_HC_chaifen.xlsx",row.names = FALSE)


#####分别导出discovery组和valuated组数据


X80_TA_dis <- X80_TA %>% filter(dis_val=="dis")
X160_HC_dis <- X160_HC %>% filter(dis_val=="dis")

X80_TA_val <- X80_TA %>% filter(dis_val=="val")
X160_HC_val <- X160_HC %>% filter(dis_val=="val")

write.xlsx(X80_TA_dis,"X80_TA_dis.xlsx",row.names = FALSE)
write.xlsx(X160_HC_dis,"X160_HC_dis.xlsx",row.names = FALSE)
write.xlsx(X80_TA_val,"X80_TA_val.xlsx",row.names = FALSE)
write.xlsx(X160_HC_val,"X160_HC_val.xlsx",row.names = FALSE)



#导入函数
source("D:/Users/wotri/Documents/R work/table/mytable_fun_liangzu_mean_sd.R")

#比较TA80与HC160
X240_total <- rbind(X80_TA[c("ID","age","sex","BMI","group", "dis_val"
)],X160_HC[c("ID","age","sex","BMI","group", "dis_val"
)])
view(X240_total)

mydata1 <- X240_total

str(mydata1)

dput(names(mydata1))

myvars1 <- c("age","sex","BMI","group")
catvars1 <- c("sex")
gvar1 <- c( "group" )

table_new1 <- mytable_fun_liangzu_mean_sd(mydata1,myvars1,catvars1,gvar1)

View(table_new1)
write.xlsx(table_new1,"x80_160.xlsx")

#比较TA中dis与val组
mydata2 <- X80_TA

str(mydata2)

dput(names(mydata2))

myvars2 <- c("age",  "BMI", "GLOB", "ESR", 
             "CRP", "sex", "NIH",  "dis_val")
catvars2 <- c("sex", "NIH")
gvar2 <- c( "dis_val")

table_new2 <- mytable_fun_liangzu_mean_sd(mydata2,myvars2,catvars2,gvar2)

View(table_new2)
write.xlsx(table_new2,"TA_dis_vs_val.xlsx")

#比较TA_dis 与 HC_dis
mydata3 <- X240_total %>% filter(dis_val=="dis")

str(mydata3)

dput(names(mydata3))

myvars3 <- c("age","sex","BMI","group")
catvars3 <- c("sex")
gvar3 <- c( "group" )

table_new3 <- mytable_fun_liangzu_mean_sd(mydata3,myvars3,catvars3,gvar3)

View(table_new3)
write.xlsx(table_new3,"TA_dis_vs_HC_dis.xlsx")


#比较TA_val 与 HC_val
mydata4 <- X240_total %>% filter(dis_val=="val")

str(mydata4)

dput(names(mydata4))

myvars4 <- c("age","sex","BMI","group")
catvars4 <- c("sex")
gvar4 <- c( "group" )

table_new4 <- mytable_fun_liangzu_mean_sd(mydata4,myvars4,catvars4,gvar4)

View(table_new4)
write.xlsx(table_new4,"TA_val_vs_HC_val.xlsx")

