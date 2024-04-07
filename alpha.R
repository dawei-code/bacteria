library(tidyverse)
library(vegan)
library(readxl)


###################
#小白鱼方法计算alpha多样性
#加载函数
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')	#Gini-Simpson 指数 
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


# plot
theme_set(theme_bw())

g <- ggplot(plot_alpha, aes(group,Richness, 
                            fill=group))
g +  geom_dotplot(binaxis='y', 
               stackdir='center', 
               dotsize = .5) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  labs(title="Richness", 
       x=NULL,
       y=NULL)


source("D:/Users/wotri/Documents/R work/table/mytable_fun_liangzu_wilcoxon.R")###两组比较函数
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

