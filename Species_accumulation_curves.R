library(vegan)
library(tidyverse)
library(xlsx)
library(readxl)

####导入数据

otu_table_raw <-as.data.frame( read_excel("otu_table_raw.xlsx"))
otu_table_raw <- column_to_rownames(otu_table_raw,var = "...1")

sample_data_raw <-as.data.frame( read_excel("sample_data_raw.xlsx"))
names(sample_data_raw)[1] <- "sample"



plot_otu <- t(otu_table_raw)
plot_otu <- rownames_to_column(as.data.frame(plot_otu),var = "sample")
plot_otu <- left_join(plot_otu,sample_data_raw[c("sample","group")],by="sample")


#计算plot_HC_data
plot_otu_HC <- plot_otu %>% filter(group=="HC") %>% select(!c(sample,group))
sp_HC <- specaccum(plot_otu_HC, method = 'random')
str(sp_HC)
ci_low_HC <- c()
ci_high_HC <- c()
for(i in 1:(nrow(sp_HC$perm))) {
  ifelse(length (unique(sp_HC$perm[i,]))==1, 
         ci_low_HC <- c(ci_low_HC,unique(sp_HC$perm[i,])), 
         ci_low_HC <- c(ci_low_HC, t.test(sp_HC$perm[i,])$conf.int[1]))
  ifelse(length (unique(sp_HC$perm[i,]))==1, 
         ci_high_HC <- c(ci_high_HC,unique(sp_HC$perm[i,])), 
         ci_high_HC <- c(ci_high_HC, t.test(sp_HC$perm[i,])$conf.int[2]))
  }
plot_HC_data <- data.frame(N_samples=sp_HC$sites,N_OTUs=sp_HC$richness,sd_OTUs=sp_HC$sd,ci_low=ci_low_HC,ci_high=ci_high_HC,Group=rep("HC",length(ci_low_HC)))

#计算plot_TA_data
plot_otu_TA <- plot_otu %>% filter(group=="TA") %>% select(!c(sample,group))
sp_TA <- specaccum(plot_otu_TA, method = 'random')
str(sp_TA)
ci_low_TA <- c()
ci_high_TA <- c()
for(i in 1:(nrow(sp_TA$perm))) {
  ifelse(length (unique(sp_TA$perm[i,]))==1, 
         ci_low_TA <- c(ci_low_TA,unique(sp_TA$perm[i,])), 
         ci_low_TA <- c(ci_low_TA, t.test(sp_TA$perm[i,])$conf.int[1]))
  ifelse(length (unique(sp_TA$perm[i,]))==1, 
         ci_high_TA <- c(ci_high_TA,unique(sp_TA$perm[i,])), 
         ci_high_TA <- c(ci_high_TA, t.test(sp_TA$perm[i,])$conf.int[2]))
}
plot_TA_data <- data.frame(N_samples=sp_TA$sites,N_OTUs=sp_TA$richness,sd_OTUs=sp_TA$sd,ci_low=ci_low_TA,ci_high=ci_high_TA,Group=rep("TA",length(ci_low_TA)))


#plot_data
plot_data <- rbind(plot_TA_data,plot_HC_data)

#作图展示
ggplot(plot_data,aes(group=Group))+
  geom_line(aes(N_samples,N_OTUs,color=Group))+
  geom_linerange(aes(x=N_samples,ymax=ci_high,ymin=ci_low,color=Group))+
  geom_ribbon(aes(x=N_samples,ymax=ci_high,ymin=ci_low,fill=Group),alpha=0.4)+
  coord_cartesian(xlim = c(0,80))


