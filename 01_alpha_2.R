#根据公司提供数据
#"reads","shannon","simpson", "chao","ace","observed_otus", 
library(tidyverse)
library(readxl)
alpha_rarefac_summary <- read_excel("alpha_rarefac.summary.xlsx")
names(alpha_rarefac_summary)[1] <- "ID"
sample_data_raw <-as.data.frame( read_excel("sample_data_raw.xlsx"))
names(sample_data_raw)[1] <- "ID"



source("../table/mytable_fun_liangzu_wilcoxon_20240405.R")###两组比较函数

mydata <- left_join(alpha_rarefac_summary,sample_data_raw,by="ID")

str(mydata)

dput(names(mydata))

myvars <- c("reads","shannon","simpson", "chao","ace","observed_otus",  "group")
catvars <- c()
gvar <- c( "group")


###############################################################################################

table_new <-mytable_fun_liangzu(mydata,myvars,catvars,gvar)

table_new$jieguo

table_new$table_fenzu_SD

# plot
#oberved otu
theme_set(theme_bw())

g <- ggplot(mydata, aes(group,observed_otus, 
                            fill=group))
g +  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Observed otu", 
       x=NULL,
       y=NULL)

#shannon
theme_set(theme_bw())

g <- ggplot(mydata, aes(group,shannon, 
                            fill=group))
g +  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Shannon", 
       x=NULL,
       y=NULL)

#Simpson
theme_set(theme_bw())

g <- ggplot(mydata, aes(group,simpson, 
                            fill=group))
g +  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Simpson", 
       x=NULL,
       y=NULL)
#Chao1

theme_set(theme_bw())

g <- ggplot(mydata, aes(group,chao, 
                            fill=group))
g +  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Chao1", 
       x=NULL,
       y=NULL)
 #ACE
theme_set(theme_bw())

g <- ggplot(mydata, aes(group,ace, 
                            fill=group))
g +  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="ACE", 
       x=NULL,
       y=NULL)
