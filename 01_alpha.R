#根据otu计算
library(tidyverse)
library(vegan)
library(readxl)


###################
#小白鱼方法计算alpha多样性
#加载函数
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]#即为observed otu
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- 1-diversity(x, index = 'simpson')	#小白鱼原始的为Gini-Simpson 指数 ，此处予以修改，改为1-Gini-Simpson
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}





otu_table_raw <-as.data.frame( read_excel("otu_table_raw.xlsx"))
otu_table_raw <- column_to_rownames(otu_table_raw,var = "...1")

sample_data_raw <-as.data.frame( read_excel("sample_data_raw.xlsx"))
names(sample_data_raw)[1] <- "ID"



plot_otu <- t(otu_table_raw)

#不包含谱系多样性，无需指定进化树；
alpha_all <- alpha(plot_otu)
alpha_all <- rownames_to_column(alpha_all,var = "ID")


#将alpha结果与sample数据集合并
plot_alpha <- left_join(alpha_all,sample_data_raw, by= "ID")

dput(names(plot_alpha))




source("../table/mytable_fun_liangzu_wilcoxon_20240405.R")###两组比较函数

mydata <- plot_alpha

str(mydata)

dput(names(mydata))

myvars <- c("Shannon","Simpson", "Chao1","ACE","Richness",  "group")
catvars <- c()
gvar <- c( "group")


###############################################################################################

table_new <-mytable_fun_liangzu(mydata,myvars,catvars,gvar)

table_new$jieguo

table_new$table_fenzu_SD

# plot
#oberved otu
theme_set(theme_bw())

g <- ggplot(plot_alpha, aes(group,Richness, 
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

g <- ggplot(plot_alpha, aes(group,Shannon, 
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

g <- ggplot(plot_alpha, aes(group,Simpson, 
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

g <- ggplot(plot_alpha, aes(group,Chao1, 
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

g <- ggplot(plot_alpha, aes(group,ACE, 
                            fill=group))
g +  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="ACE", 
       x=NULL,
       y=NULL)