#本文需要用到的 R 包
library(vegan)	#用于计算 Shannon 熵指数、Simpson 指数、Chao1 指数、ACE 指数等，同时用于抽样
library(picante)	#用于计算 PD_whole_tree，若不计算它就无需加载。事实上，picante 包加载时默认同时加载 vegan
library(tidyverse)	#用于 ggplot2 作图
library(doBy)	#用于分组统计（第 97 行使用到）
library(ggalt)	#用于绘制拟合曲线（第 138 行使用到）
library(BiodiversityR)	#用于绘制 Rank-abundance 曲线（第 160 行使用到）
library(xlsx)
library(readxl)

##############################
##导入自定义函数
#计算多种 Alpha 多样性指数，结果返回至向量

alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)	#丰富度指数
  else if (method == 'chao1') result <- estimateR(x)[2, ]	#Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]	#ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)	#Shannon 指数
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')	#Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)	#Pielou 均匀度
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)	#goods_coverage
  else if (method == 'pd' & !is.null(tree)) {	#PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}

#根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表,rare参数指终点
alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
  x_nrow <- nrow(x)
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step)
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
    
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  
  names(alpha_rare) <- rownames(x)
  alpha_rare
}



####导入数据

otu_table_raw <-as.data.frame( read_excel("otu_table_raw.xlsx"))
otu_table_raw <- column_to_rownames(otu_table_raw,var = "...1")

sample_data_raw <-as.data.frame( read_excel("sample_data_raw.xlsx"))
names(sample_data_raw)[1] <- "sample"



plot_otu <- t(otu_table_raw)

###平时常说的稀释曲线（rarefaction curve），绝大多数情况下等同于Richness指数（物种丰富度指数）的Alpha多样性曲线
richness_curves <- alpha_curves(plot_otu, step =1000, method = 'richness')
#获得 ggplot2 作图文件
plot_richness <- data.frame()
for (i in names(richness_curves)) {
  richness_curves_i <- (richness_curves[[i]])
  richness_curves_i <- data.frame(rare = names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
  plot_richness <- rbind(plot_richness, richness_curves_i)
}

rownames(plot_richness) <- NULL
plot_richness$rare <- as.numeric(plot_richness$rare)
plot_richness$alpha <- as.numeric(plot_richness$alpha)

plot_richness <- left_join(plot_richness, sample_data_raw,by="sample")

#ggplot2 作图
ggplot(plot_richness, aes(rare, alpha, group = sample,color=group)) +
  
  labs(x = 'Number of sequences', y = 'Richness', color = NULL) +
  geom_line()+
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = min(rowSums(plot_otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 100000, 10000), labels = as.character(seq(0, 100000, 10000)))

##多计算几次以获取均值 ± 标准差，然后再展示出也是一个不错的选择
#重复抽样 5 次
plot_richness <- data.frame()

for (n in 1:5) {
  richness_curves <- alpha_curves(plot_otu, step = 1000, method = 'richness')
  
  for (i in names(richness_curves)) {
    richness_curves_i <- (richness_curves[[i]])
    richness_curves_i <- data.frame(rare = names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
    plot_richness <- rbind(plot_richness, richness_curves_i)
  }
}

#计算均值 ± 标准差（doBy 包中的 summaryBy() 函数）
plot_richness_stat <- summaryBy(alpha~sample+rare, plot_richness, FUN = c(mean, sd))
plot_richness_stat$rare <- as.numeric(plot_richness_stat$rare)
plot_richness_stat[which(plot_richness_stat$rare == 0),'alpha.sd'] <- NA

plot_richness_stat <- plot_richness_stat[order(plot_richness_stat$sample, plot_richness_stat$rare), ]
plot_richness_stat <- left_join(plot_richness_stat, sample_data_raw,by="sample")
#ggplot2 作图
ggplot(plot_richness_stat, aes(rare, alpha.mean, group=sample,color=group)) +
  geom_xspline(spline_shape=1) +
  labs(x = 'Number of sequences', y = 'Richness', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = min(rowSums(plot_otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 100000, 10000), labels = as.character(seq(0, 100000, 10000)))

##对于 Shannon 指数等，方法类似
#以 2000 步长（step=2000）为例统计每个稀释梯度下的 Shannon 指数，Shannon 公式的对数底数默认为 e，若有需要可更改（例如 2）
shannon_curves <- alpha_curves(plot_otu, step = 2000, method = 'shannon',base = exp(1))

#获得 ggplot2 作图文件（略，参见上述）
#ggplot2 作图（略，参见上述）

##若简单的“geom_line()”样式波动幅度过大，不平滑等，可以尝试拟合曲线的样式
#获得作图数据。前面多生成一个点，使得 Shannon 拟合曲线更加平滑（你把 shannon_curves1 注释掉就知道我说的啥了）
shannon_curves1 <- alpha_curves(plot_otu, step = 200, rare = 200, method = 'shannon')
shannon_curves2 <- alpha_curves(plot_otu, step = 2000, method = 'shannon')
shannon_curves <- c(shannon_curves1, shannon_curves2)

plot_shannon <- data.frame()
for (i in 1:length(shannon_curves)) {
  shannon_curves_i <- shannon_curves[[i]]
  shannon_curves_i <- data.frame(rare = names(shannon_curves_i), alpha = shannon_curves_i, sample = names(shannon_curves)[i], stringsAsFactors = FALSE)
  plot_shannon <- rbind(plot_shannon, shannon_curves_i)
}

rownames(plot_shannon) <- NULL
plot_shannon$rare <- as.numeric(plot_shannon$rare)
plot_shannon$alpha <- as.numeric(plot_shannon$alpha)
plot_shannon <- left_join(plot_shannon, sample_data_raw,by="sample")

plot_shannon <- plot_shannon[order(plot_shannon$sample, plot_shannon$rare), ]

#ggplot2 作图（使用到 ggalt 包的 geom_xspline() 绘制平滑拟合线）
ggplot(plot_shannon, aes(rare, alpha, group=sample,color = group)) +
  geom_smooth(se=FALSE) +
  labs(x = 'Number of sequences', y = 'Shannon', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = min(rowSums(plot_otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 30000, 5000), labels = as.character(seq(0, 30000, 5000)))



##Rank-abundance 曲线
#统计（BiodiversityR 包 rankabundance() 实现 OTU 排序）
otu_relative <- plot_otu / rowSums(plot_otu)
rank_dat <- data.frame()
for (i in rownames(otu_relative)) {
  rank_dat_i <- data.frame(rankabundance(subset(otu_relative, rownames(otu_relative) == i), digits = 6))[1:2]
  rank_dat_i$sample <- i
  rank_dat <- rbind(rank_dat, rank_dat_i)
}
rank_dat <- subset(rank_dat, rank_dat$abundance != 0)

rank_dat <- left_join(rank_dat,sample_data_raw,by="sample")

#ggplot2 作图
ggplot(rank_dat, aes(rank, log(abundance, 10), group = sample,color=group)) +
  geom_line() +
  labs(x = 'OTUs rank', y = 'Relative adundance (%)', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  scale_y_continuous(breaks = 0:-5, labels = c('100', '10', '1', '0.1', '0.01', '0.001'), limits = c(-5, 0))

